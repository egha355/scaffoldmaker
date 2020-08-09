"""
Generates a 3-D unit cylinder mesh with variable numbers of elements around, along and
through wall.
"""

from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.torso_terms import get_torso_term
from scaffoldmaker.utils import geometry
from scaffoldmaker.utils.shieldmesh import ShieldMesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.interpolation import computeCubicHermiteDerivativeScaling, getCubicHermiteArcLength, interpolateSampleCubicHermite, \
    sampleCubicHermiteCurves, sampleCubicHermiteCurvesSmooth, smoothCubicHermiteDerivativesLine
from opencmiss.utils.zinc.finiteelement import getMaximumElementIdentifier, getMaximumNodeIdentifier


class MeshType_3d_solidcylinder9(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Solid Cylinder 9'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Major radius' : 2.0,
            'Minor radius' : 1.0,
            'Height' : 5.0,
            'Number of elements across torso' : 8,
            'Number of elements up torso' : 5,
            'Number of elements along torso' : 1,
            'Number of elements for rim' : 0,
            'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Major radius',
            'Minor radius',
            'Height',
            'Number of elements across torso',
            'Number of elements up torso',
            'Number of elements along torso',
            'Number of elements for rim',
            'Use cross derivatives',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements across torso',
            'Number of elements along torso']:
            if options[key] < 1:
                options[key] = 1
        if (options['Number of elements across torso'] < 2) :
            options['Number of elements across torso'] = 8

    @staticmethod
    def createEllipsePoints(centre, majorAxis, minorAxis, elementsCountAround, height):
        '''
        Generate a set of points and derivatives for an ellipse
        starting at pole majorAxis from centre.
        :param centre: Centre of full ellipse.
        :param majorAxis: Vector in direction of starting major radius, magnitude is ellipse major radius.
        :param minorAxis: Vector normal to major axis, magnitude is ellipse minor axis length.
        :param height: Height of arc of ellipsoid from starting point along majorAxis.
        :return: Lists nx, nd1. Ordered fastest around, starting at major radius.
        '''
        nx = []
        nd1 = []
        magMajorAxis = vector.magnitude(majorAxis)
        magMinorAxis = vector.magnitude(minorAxis)
        unitMajorAxis = vector.normalise(majorAxis)
        unitMinorAxis = vector.normalise(minorAxis)
        useHeight = min(max(0.0, height), 2.0 * magMajorAxis)
        totalRadians = geometry.getEllipseRadiansToX(magMajorAxis, 0.0, magMajorAxis - useHeight,
                                              initialTheta=0.5 * math.pi * useHeight / magMajorAxis)
        radians = 0.0
        arcLengthUp = geometry.getEllipseArcLength(magMajorAxis, magMinorAxis, radians, totalRadians)
        elementsCountUp = elementsCountAround//2
        elementArcLength = arcLengthUp / elementsCountUp
        # radiansPerElementAround = 2.0 * math.pi / elementsCountAround
        for n1 in range(elementsCountUp+1):
            for n2 in [-1,1]:
                cosRadians = math.cos(radians)
                sinRadians = n2*math.sin(radians)
                nx.append(
                    [ (centre[c] + cosRadians * majorAxis[c] +sinRadians * minorAxis[c]) for c in range(3) ])

                ndab = vector.setMagnitude([-sinRadians * magMajorAxis, cosRadians * magMinorAxis], elementArcLength)
                nd1.append(
                        [(ndab[0] * unitMajorAxis[c] + ndab [1] * unitMinorAxis[c]) for c in range(3)])
            radians = geometry.updateEllipseAngleByArcLength(magMajorAxis, magMinorAxis, radians, elementArcLength)
        return nx, nd1

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        majorRadius = options['Major radius']
        minorRadius = options['Minor radius']
        height = options['Height']
        elementsCountAcross = options['Number of elements across torso']
        elementsCountUp = options['Number of elements up torso']
        elementsCountRim = options['Number of elements for rim']
        elementsCountAlong = options['Number of elements along torso']
        elementsCountUpRegular = elementsCountUp - 2 - elementsCountRim
        elementsCountAcrossNonRim = elementsCountAcross - 2*elementsCountRim
        elementsCountAround = 2 * elementsCountUpRegular + elementsCountAcrossNonRim
        useCrossDerivatives = options['Use cross derivatives']


        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)


        btGroup = AnnotationGroup(region, get_torso_term("anterior of torso"))
        annotationGroups = [ btGroup ]

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        mesh = fm.findMeshByDimension(3)
        cache = fm.createFieldcache()

        nx, nd1 = cls.createEllipsePoints(
            [0.0, 0.0, 0.0], [0.0 , 0.0, -majorRadius], [0.0, minorRadius, 0.0], elementsCountAround, majorRadius)
        normalToEllipse = vector.normalise(vector.crossproduct3([0.0 , 0.0, -majorRadius], [0.0, minorRadius, 0.0]))



        # The bottom curve node coordinates and derivatives
        tShield = ShieldMesh( elementsCountAcross, elementsCountUp, elementsCountRim, None, elementsCountAlong)

        tix = []
        tox = []
        tnx = []
        tid1 = []
        tod1 = []
        tid2 = []
        tod2 = []
        tid3 = []
        tod3 = []
        gn = elementsCountAround
        for n in range(elementsCountAround+1):
            tix.append(nx[gn])
            tid2.append(nd1[gn])
            tid3.append([height*normalToEllipse[c] for c in range(3)])
            left, right = (n<elementsCountAround//2), (n>=elementsCountAround//2)
            gn += 2*(right-left)+(gn==0)

            tod3.append(tid3[n])
            tid1.append(
                vector.normalise(vector.crossproduct3(tid2[n], tid3[n]))
            )
            tod2.append(tid2[n])
            tod1.append(tid1[n])

        for n3 in range(1,tShield.elementsCountThroughWall+1):
            tox = []
            for n in range(elementsCountAround + 1):
                tox.append(
                    [ (tix[n][c] + n3*height * normalToEllipse[c]) for c in range(3) ])
            tnx.append(tox)

        btx  = tShield.px
        btd1 = tShield.pd1
        btd2 = tShield.pd2
        btd3 = tShield.pd3
        for n in range(tShield.elementsCountAroundFull + 1):
            n1, n2 = tShield.convertRimIndex(n)
            for n3 in range(tShield.elementsCountThroughWall+1):
                if n3 == 0:
                    tx, td1, td2, td3 = tix, tid1, tid2, tid3
                else:
                    tx, td1, td2, td3 = tnx[n3-1], tod1, tod2, tod3
                btx[n3][n2][n1] = tx[n]
                if n2 > tShield.elementsCountRim:  # regular rows
                    if n1 < tShield.elementsCountAcross:
                        btd1[n3][n2][n1] = [ -d for d in td1[n] ]
                        btd2[n3][n2][n1] = [ -d for d in td2[n] ]
                    else:
                        btd1[n3][n2][n1] = td1[n]
                        btd2[n3][n2][n1] = td2[n]
                else:  # around rim
                    btd1[n3][n2][n1] = td2[n]
                    btd2[n3][n2][n1] = [ -d for d in td1[n] ]
                btd3[n3][n2][n1] = td3[n]

        rcx = []
        tmdx = btx[0][0][(tShield.elementsCountAcross)//2]
        tmdd2 = btd2[0][0][(tShield.elementsCountAcross)//2]
        tmux = [0.5*(btx[0][tShield.elementsCountUp][0][c]+btx[0][tShield.elementsCountUp][tShield.elementsCountAcross][c]) for c in range(3)]
        rcx.append(tmdx)
        rcx.append(tmux)
        rcd2 = [tmdd2,tmdd2]
        rscx, rscd2 = sampleCubicHermiteCurves(rcx, rcd2, tShield.elementsCountUp, lengthFractionStart=0.667, arcLengthDerivatives = True)[0:2]

        # get d1, d3
        rscd1 = []
        rscd3 = []
        for n in range(len(rscx)):
            d1 = vector.normalise([tix[tShield.elementsCountAroundFull][c]-tix[0][c] for c in range(3)])
            d3 = vector.normalise(vector.crossproduct3(d1, rscd2[n]))
            rscd1.append(d1)
            rscd3.append(d3)

        # across regular rows of Shield: get d1, initial d2
        for n2 in range(tShield.elementsCountRim + 2, tShield.elementsCountUp + 1):
            btx[0][n2], btd1[0][n2], pe, pxi, psf = sampleCubicHermiteCurves(
                [ btx[0][n2][0], rscx[n2], btx[0][n2][-1] ], [ btd1[0][n2][0], rscd1[n2], btd1[0][n2][-1] ], tShield.elementsCountAcross,
                lengthFractionStart=0.667, lengthFractionEnd=0.667, arcLengthDerivatives = True)
            btd2[0][n2] = interpolateSampleCubicHermite([btd2[0][n2][0], rscd2[n2], btd2[0][n2][-1]], [[0.0, 0.0, 0.0]] * 3, pe, pxi, psf)[0]

        # up regular columns of shield: get d2, initial d1 below regular rows
        for n1 in range(2, tShield.elementsCountAcross - 1):
            tx, td2, pe, pxi, psf = sampleCubicHermiteCurves(
                [ btx[0][0][n1], btx[0][2][n1] ], [ btd2[0][0][n1], btd2[0][2][n1] ], 2, lengthFractionStart=0.667, arcLengthDerivatives = True)  # GRC fudge factor rvSulcusEdgeFactor
            for n2 in range(3, tShield.elementsCountUp + 1):
                tx .append(btx [0][n2][n1])
                td2.append(btd2[0][n2][n1])
            td2 = smoothCubicHermiteDerivativesLine(tx, td2, fixStartDirection = True)
            td1 = interpolateSampleCubicHermite([ btd1[0][0][n1], btd1[0][2][n1] ], [ [ 0.0, 0.0, 0.0 ] ]*2, pe, pxi, psf)[0]
            for n2 in range(tShield.elementsCountUp + 1):
                if n2 < 2:
                    btx [0][n2][n1] = tx [n2]
                    btd1[0][n2][n1] = td1[n2]
                btd2[0][n2][n1] = td2[n2]

        tShield.getTriplePoints(n3=0)
        n1b = 1
        m1a = tShield.elementsCountAcross
        m1b = m1a - 1
        m1c = m1a - 2
        n2b = 1

        # smooth RV freewall row 1
        btd1[0][n2b][n1b:m1a] = smoothCubicHermiteDerivativesLine(btx[0][n2b][n1b:m1a], btd1[0][n2b][n1b:m1a])

        # smooth RV columns 1, -2
        for n1 in [ 1, -2 ]:
            tx = []
            td2 = []
            for n2 in range(1, tShield.elementsCountUp + 1):
                tx .append(btx [0][n2][n1])
                td2.append(btd2[0][n2][n1])
            td2 = smoothCubicHermiteDerivativesLine(tx, td2)
            for n in range(tShield.elementsCountUp):
                btd2[0][n + 1][n1] = td2[n]

        tShield.smoothDerivativesToTriplePoints(n3=0, fixAllDirections=True)

        # get outer d3 and inner x, d3
        for n2 in range(tShield.elementsCountUp + 1):
            for n1 in range(tShield.elementsCountAcross + 1):
                if btd1[0][n2][n1]:
                    btd3[0][n2][n1] = btd3[1][n2][n1] = vector.setMagnitude(vector.crossproduct3(btd1[0][n2][n1], btd2[0][n2][n1]), height)
                    btx [1][n2][n1] = [ (btx [0][n2][n1][c] + btd3[0][n2][n1][c]) for c in range(3) ]

        # get inner d1, d2
        # row 1
        btd1[1][n2b][n1b:m1a] = smoothCubicHermiteDerivativesLine(btx[1][n2b][n1b:m1a], btd1[0][n2b][n1b:m1a], fixAllDirections = True)
        # regular rows 2+
        for n2 in range(2, tShield.elementsCountUp + 1):
            btd1[1][n2] = smoothCubicHermiteDerivativesLine(btx[1][n2], btd1[0][n2], fixAllDirections = True)
        # columns
        for n1 in range(n1b, m1a):
            startn2 = 1 if (n1 in [n1b, m1b]) else 0
            tx  = []
            td2 = []
            for n2 in range(startn2, tShield.elementsCountUp + 1):
                tx .append(btx [1][n2][n1])
                td2.append(btd2[0][n2][n1])
            td2 = smoothCubicHermiteDerivativesLine(tx, td2, fixAllDirections = True)
            for n2 in range(startn2, tShield.elementsCountUp + 1):
                btd2[1][n2][n1] = td2[n2 - startn2]

        # fix inner derivatives leading to triple points
        # first copy d2 from outer to inner
        for n1 in [ n1b, m1b ]:
            for n2 in range(0, n2b):
                btd2[1][n2][n1] = btd2[0][n2][n1]
        tShield.smoothDerivativesToTriplePoints(n3=1, fixAllDirections=True)

        btd2[1][0][1] = smoothCubicHermiteDerivativesLine([ btx[1][0][1], btx[1][1][1] ], [ btd2[1][0][1], [ (btd1[1][1][1][c] + btd2[1][1][1][c]) for c in range(3) ] ],
                                                          fixEndDerivative = True, fixStartDirection = True)[0]
        btd2[1][0][m1b] = smoothCubicHermiteDerivativesLine([ btx[1][0][m1b], btx[1][1][m1b] ], [ btd2[1][0][m1b], [ (-btd1[1][1][m1b][c] + btd2[1][1][m1b][c]) for c in range(3) ] ],
                                                          fixEndDerivative = True, fixStartDirection = True)[0]


        temx = []
        if tShield.elementsCountThroughWall>1:
            for n2 in range(tShield.elementsCountUp + 1):
                for n3 in range(2, tShield.elementsCountThroughWall + 1):
                    for n1 in range(tShield.elementsCountAcross + 1):
                        if tShield.px[0][n2][n1]:
                            temx = [tShield.px[0][n2][n1][c] + n3*height*normalToEllipse[c] for c in range(3)]
                            tShield.px[n3][n2][n1]=temx
                            tShield.pd1[n3][n2][n1]=tShield.pd1[0][n2][n1]
                            tShield.pd2[n3][n2][n1]=tShield.pd2[0][n2][n1]
                            tShield.pd3[n3][n2][n1]=tShield.pd3[0][n2][n1]

        #################
        # Create nodes
        #################

        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        nodeIdentifier = startNodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        nodeIdentifier = tShield.generateNodes(fm, coordinates, nodeIdentifier)

        #################
        # Create elements
        #################

        btMeshGroup = btGroup.getMeshGroup(mesh)
        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)
        elementIdentifier = tShield.generateElements(fm, coordinates, elementIdentifier, [btMeshGroup])

        fm.endChange()
        return annotationGroups

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        """
        if not options['Refine']:
            cls.generateBaseMesh(region, options)
            return

        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        baseRegion = region.createRegion()
        cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong, refineElementsCountThroughWall)

    @classmethod
    def defineFaceAnnotations(cls, region, options, annotationGroups):
        """
        Add face annotation groups from the highest dimension mesh.
        Must have defined faces and added subelements for highest dimension groups.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for top-level elements.
        New face annotation groups are appended to this list.
        """
        # create 2d surface mesh groups
        fm = region.getFieldmodule()
        anteriorGroup = getAnnotationGroupForTerm(annotationGroups, get_torso_term("anterior of torso"))
        posteriorGroup = getAnnotationGroupForTerm(annotationGroups, get_torso_term("posterior of torso"))

        mesh2d = fm.findMeshByDimension(2)
        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))

        is_anterior = anteriorGroup.getFieldElementGroup(mesh2d)
        is_anterior_faces = fm.createFieldAnd(is_anterior, is_exterior_face_xi3_1)


        is_posterior = posteriorGroup.getFieldElementGroup(mesh2d)
        is_posterior_faces = fm.createFieldAnd(is_posterior, is_exterior_face_xi3_1)

        anteriorFaces = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_torso_term("anterior of torso"))
        anteriorFaces.getMeshGroup(mesh2d).addElementsConditional(is_anterior_faces)
        posteriorFaces = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_torso_term("posterior of torso"))
        posteriorFaces.getMeshGroup(mesh2d).addElementsConditional(is_posterior_faces)
