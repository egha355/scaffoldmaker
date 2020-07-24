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
# from scaffoldmaker.utils.sheildmesh import ShieldMesh


class MeshType_3d_solidcylinder7(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Solid Cylinder 7'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Major radius' : 263,
            'Minor radius' : 190,
            'Height' : 600,
            'Number of elements around' : 16,
            'Number of elements along' : 1,
            'Number of elements through wall' : 1,
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
            'Number of elements around',
            'Number of elements along',
            'Number of elements through wall',
            'Use cross derivatives',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements along',
            'Number of elements through wall',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if (options['Number of elements around'] < 2) :
            options['Number of elements around'] = 2


    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        majorRadius = options['Major radius']
        minorRadius = options['Minor radius']
        height = options['Height']
        elementsCountAround = options['Number of elements around']
        elementsCountAlong = options['Number of elements along']
        elementsCountThroughWall = options['Number of elements through wall']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        # myShield = ShieldMesh(4, 2, 0, None)
        # # myShield.getTriplePoints(5)
        # myShield.generateNodes(fm, coordinates, 1)
        # myShield.generateElements(fm, coordinates, 1, [])

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        mesh = fm.findMeshByDimension(3)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        eft2 = tricubichermite.createEftWedgeRadial(0 * 100, 1* 100)
        elementtemplate2 = mesh.createElementtemplate()
        elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate2.defineField(coordinates, -1, eft2)




        anteriorGroup = AnnotationGroup(region, get_torso_term("anterior of torso"))
        posteriorGroup = AnnotationGroup(region, get_torso_term("posterior of torso"))
        annotationGroups = [anteriorGroup, posteriorGroup]

        anteriorMeshGroup = anteriorGroup.getMeshGroup(mesh)
        posteriorMeshGroup = posteriorGroup.getMeshGroup(mesh)






        cache = fm.createFieldcache()

        # create nodes
        nodeIdentifier = 1
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0 , 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, height / elementsCountAlong ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]

        radiansPerElementAround = 2 * math.pi / elementsCountAround
        elementsCountOnLine = int(elementsCountAround/4)

        for nz in range(elementsCountAlong+1):
            x[2] = nz * height/elementsCountAlong+800
            for nr in range(elementsCountThroughWall):
                for nt in range(elementsCountAround):
                    radiansAround = nt * radiansPerElementAround
                    # radiansAround = geometry.updateEllipseAngleByArcLength(majorRadius, minorRadius, 0, nt*arcLengthPerElementAround)
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    x[0] = majorRadius*cosRadiansAround
                    x[1] = minorRadius*sinRadiansAround
                    # DTS = 1.0/math.sqrt(majorRadius/minorRadius*x[1]*majorRadius/minorRadius*x[1]+minorRadius/majorRadius*x[0]*minorRadius/majorRadius*x[0])
                    # dx_ds1 = [-majorRadius*sinRadiansAround*DTS*arcLengthPerElementAround, minorRadius*cosRadiansAround*DTS*arcLengthPerElementAround, 0.0]
                    dx_ds1 = [-majorRadius*sinRadiansAround*radiansPerElementAround, minorRadius*cosRadiansAround*radiansPerElementAround, 0.0]

                    dx_ds3 = [(nt%int(elementsCountAround/2) == 0)*majorRadius*(1-math.cos(radiansPerElementAround)), x[1], 0.0]
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    nodeIdentifier += 1


        # line nodes
        for nz in range(elementsCountAlong+1):
            for nt in range(2*elementsCountOnLine-1):
                radiansAround = (nt+1) * radiansPerElementAround
                radiansAroundNext = radiansAround - radiansPerElementAround
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)

                x = [majorRadius*cosRadiansAround, 0.0,  nz * height / elementsCountAlong+800]
                dx_ds1 = [majorRadius*(math.cos(radiansAroundNext)-cosRadiansAround), 0.0, 0.0]
                dx_ds3 = [0.0, -minorRadius*sinRadiansAround, 0.0]
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1



        # create cubic elements

        elementIdentifier = 1
        nowL = 2*elementsCountOnLine-1
        nowT = (elementsCountAlong+1)*elementsCountAround
        for ez in range(elementsCountAlong+1):
            for et in range(1, nowL):
                bni6 = (et+1)+ez*elementsCountAround
                bni2 = nowT + et + ez * nowL
                nodeIdentifiers=[bni2+1,bni2,bni2+1+nowL,bni2+nowL,bni6+1,bni6,bni6+1+elementsCountAround,bni6+elementsCountAround]
                element = mesh.createElement(elementIdentifier, elementtemplate)
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                elementIdentifier += 1
                posteriorMeshGroup.addElement(element)
            for et in range(1, nowL):
                bni6 = (ez+1)*elementsCountAround-et+1
                bni2 = nowT + et + ez * nowL
                nodeIdentifiers=[bni2+1,bni2,bni2+1+nowL,bni2+nowL,
                                 bni6-1,bni6,bni6-1+elementsCountAround,bni6+elementsCountAround]

                element = mesh.createElement(elementIdentifier, elementtemplate)
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                elementIdentifier += 1
                anteriorMeshGroup.addElement(element)


        # create wedge elements
        for ez in range(elementsCountAlong+1):
            bni3 = ez * elementsCountAround + 1
            bni1 = nowT + ez * nowL + 1
            nodeIdentifiers2 = [bni1, bni1+nowL, bni3, bni3+1, bni3+elementsCountAround, bni3+1+elementsCountAround]
            element = mesh.createElement(elementIdentifier, elementtemplate2)
            result = element.setNodesByIdentifier(eft2, nodeIdentifiers2)
            elementIdentifier += 1

            bni3 = ez * elementsCountAround + int(elementsCountAround/2)
            bni1 = nowT + ez * nowL + int(elementsCountAround/2)-1
            nodeIdentifiers2 = [bni1, bni1+nowL, bni3, bni3+1, bni3+elementsCountAround, bni3+1+elementsCountAround]
            element = mesh.createElement(elementIdentifier, elementtemplate2)
            result = element.setNodesByIdentifier(eft2, nodeIdentifiers2)
            elementIdentifier += 1

            bni3 = (ez+1) * elementsCountAround
            bni1 = nowT + ez * nowL + 1
            nodeIdentifiers2 = [bni1, bni1+nowL, bni3, bni3-elementsCountAround+1, bni3+elementsCountAround, bni3+1]
            element = mesh.createElement(elementIdentifier, elementtemplate2)
            result = element.setNodesByIdentifier(eft2, nodeIdentifiers2)
            elementIdentifier += 1

            bni3 = ez * elementsCountAround + nowL+2
            bni1 = nowT + ez * nowL + nowL
            nodeIdentifiers2 = [bni1, bni1+nowL, bni3, bni3+1, bni3+elementsCountAround, bni3+1+elementsCountAround]
            element = mesh.createElement(elementIdentifier, elementtemplate2)
            result = element.setNodesByIdentifier(eft2, nodeIdentifiers2)
            elementIdentifier += 1


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
