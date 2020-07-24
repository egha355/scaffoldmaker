"""
Generates 3-D mesh of left and right ventricles below base plane.
Variant using collapsed/wedge elements at septum junction.
"""

from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldNodeGroup, findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.heart_terms import get_heart_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds
from scaffoldmaker.utils.geometry import getApproximateEllipsePerimeter, getEllipseArcLength, updateEllipseAngleByArcLength
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class MeshType_3d_heartventricles1(Scaffold_base):
    '''
    Generates 3-D mesh of left and right ventricles below base plane.
    '''

    @staticmethod
    def getName():
        return '3D Heart Ventricles 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Mouse 1',
            'Pig 1',
            'Rat 1',
            'Unit Human 1',
            'Unit Mouse 1',
            'Unit Pig 1',
            'Unit Rat 1']

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        isHuman = 'Human' in parameterSetName
        isMouse = 'Mouse' in parameterSetName
        isPig = 'Pig' in parameterSetName
        isRat = 'Rat' in parameterSetName
        notUnitScale = 'Unit' not in parameterSetName
        options = {}
        options['Number of elements around LV free wall'] = 5
        options['Number of elements around RV free wall'] = 7
        options['Number of elements up LV apex'] = 1
        options['Number of elements up RV'] = 3
        options['Unit scale'] = 1.0
        options['Interventricular sulcus derivative factor'] = 0.5
        options['LV outer height'] = 1.0
        options['LV outer diameter'] = 0.95
        options['LV free wall thickness'] = 0.14
        options['LV apex thickness'] = 0.08
        options['RV inner height fraction'] = 0.85
        options['RV arc around degrees'] = 145.0
        options['RV arc apex fraction'] = 0.5
        options['RV free wall thickness'] = 0.05
        options['RV width'] = 0.3
        options['RV width growth factor'] = 0.7
        options['RV side extension'] = 0.12
        options['RV side extension growth factor'] = 0.5
        options['Ventricular septum thickness'] = 0.12
        options['Ventricular septum base radial displacement'] = 0.1
        options['Use cross derivatives'] = False
        options['Refine'] = False
        options['Refine number of elements surface'] = 4
        options['Refine number of elements through LV wall'] = 1
        options['Refine number of elements through wall'] = 1
        if isHuman:
            if 'Unit' not in parameterSetName:
                options['Unit scale'] = 80.0
        elif isMouse or isRat:
            if notUnitScale:
                options['Unit scale'] = 5.0 if isMouse else 12.0
            options['Interventricular sulcus derivative factor'] = 0.8
            options['LV outer height'] = 0.85
            options['LV outer diameter'] = 0.85
            options['LV free wall thickness'] = 0.12
            options['LV apex thickness'] = 0.05
            options['RV inner height fraction'] = 0.9
            options['RV arc around degrees'] = 140.0
            options['RV arc apex fraction'] = 0.5
            options['RV free wall thickness'] = 0.05
            options['RV width'] = 0.28
            options['RV width growth factor'] = 0.8
            options['RV side extension'] = 0.05
            options['RV side extension growth factor'] = 0.1
            options['Ventricular septum thickness'] = 0.15
            options['Ventricular septum base radial displacement'] = 0.0
        elif isPig:
            options['Number of elements up LV apex'] = 1
            options['Number of elements up RV'] = 3
            if 'Unit' not in parameterSetName:
                options['Unit scale'] = 80.0
            options['LV outer height'] = 0.9
            options['LV free wall thickness'] = 0.17
            options['LV apex thickness'] = 0.07
            options['RV inner height fraction'] = 0.65
            options['RV arc around degrees'] = 140.0
            options['RV free wall thickness'] = 0.06
            options['RV width growth factor'] = 0.4
            options['RV side extension'] = 0.1
            options['RV side extension growth factor'] = 0.4
            options['Ventricular septum thickness'] = 0.13
        elif 'Rat' in parameterSetName:
            if 'Unit' not in parameterSetName:
                options['Unit scale'] = 12.0
            options['LV outer height'] = 0.9
            options['LV apex thickness'] = 0.08
            options['RV width'] = 0.35
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around LV free wall',
            'Number of elements around RV free wall',
            'Number of elements up LV apex',
            'Number of elements up RV',
            'Unit scale',
            'Interventricular sulcus derivative factor',
            'LV outer height',
            'LV outer diameter',
            'LV free wall thickness',
            'LV apex thickness',
            'RV inner height fraction',
            'RV arc around degrees',
            'RV arc apex fraction',
            'RV free wall thickness',
            'RV width',
            'RV width growth factor',
            'RV side extension',
            'RV side extension growth factor',
            'Ventricular septum thickness',
            'Ventricular septum base radial displacement',
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through LV wall',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        '''
        :return:  True if dependent options changed, otherwise False.
        '''
        dependentChanges = False
        for key in [
            'Refine number of elements surface',
            'Refine number of elements through LV wall',
            'Refine number of elements through wall',
            'Number of elements up LV apex']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Number of elements around LV free wall',
            'Number of elements up RV']:
            if options[key] < 2:
                options[key] = 2
        for key in [
            'Number of elements around RV free wall']:
            if options[key] < 3:
                options[key] = 3
        for key in [
            'Unit scale',
            'LV outer height',
            'LV outer diameter',
            'LV free wall thickness',
            'LV apex thickness',
            'RV free wall thickness',
            'RV width',
            'RV side extension',
            'Ventricular septum thickness']:
            if options[key] < 0.0:
                options[key] = 0.0
        for key in [
            'RV inner height fraction',
            'RV width growth factor',
            'RV side extension growth factor']:
            if options[key] < 0.001:
                options[key] = 0.001
            elif options[key] > 0.999:
                options[key] = 0.999
        if options['RV arc around degrees'] < 45.0:
            options['RV arc around degrees'] = 45.0
        elif options['RV arc around degrees'] > 270.0:
            options['RV arc around degrees'] = 270.0
        for key in [
            'Interventricular sulcus derivative factor',
            'RV arc apex fraction']:
            if options[key] < 0.1:
                options[key] = 0.1
            elif options[key] > 1.0:
                options[key] = 1.0
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        elementsCountAroundLVFreeWall = options['Number of elements around LV free wall']
        elementsCountAroundRVFreeWall = options['Number of elements around RV free wall']
        pointsCountAroundOuter = elementsCountAroundLV = elementsCountAroundLVFreeWall + elementsCountAroundRVFreeWall
        elementsCountAroundVSeptum = elementsCountAroundRVFreeWall
        elementsCountUpLVApex = options['Number of elements up LV apex']
        elementsCountUpRV = options['Number of elements up RV']
        elementsCountUpLV = elementsCountUpLVApex + elementsCountUpRV
        unitScale = options['Unit scale']
        ivSulcusDerivativeFactor = options['Interventricular sulcus derivative factor']

        lvOuterHeight = unitScale*options['LV outer height']
        lvOuterRadius = unitScale*0.5*options['LV outer diameter']
        lvFreeWallThickness = unitScale*options['LV free wall thickness']
        lvApexThickness = unitScale*options['LV apex thickness']
        lvInnerHeight = lvOuterHeight - lvApexThickness
        lvInnerRadius = lvOuterRadius - lvFreeWallThickness
        rvInnerHeightFraction = options['RV inner height fraction']
        rvArcAroundBaseRadians = math.radians(options['RV arc around degrees'])
        rvArcApexFraction = options['RV arc apex fraction']
        rvFreeWallThickness = unitScale*options['RV free wall thickness']
        rvWidth = unitScale*options['RV width']
        rvWidthGrowthFactor = options['RV width growth factor']
        rvSideExtension = unitScale*options['RV side extension']
        rvSideExtensionGrowthFactor = options['RV side extension growth factor']
        vSeptumThickness = unitScale*options['Ventricular septum thickness']
        vSeptumBaseRadialDisplacement = unitScale*options['Ventricular septum base radial displacement']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)
        cache = fm.createFieldcache()

        mesh = fm.findMeshByDimension(3)

        lvGroup = AnnotationGroup(region, get_heart_term("left ventricle myocardium"))
        rvGroup = AnnotationGroup(region, get_heart_term("right ventricle myocardium"))
        vSeptumGroup = AnnotationGroup(region, get_heart_term("interventricular septum"))
        annotationGroups = [ lvGroup, rvGroup, vSeptumGroup ]

        # annotation fiducial points
        markerGroup = findOrCreateFieldGroup(fm, "marker")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)

        #################
        # Create nodes
        #################
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        nodetemplateApex = nodetemplate

        nodeIdentifier = 1

        # LV nodes

        # get distance from outside of round outer LV to outside of RV at base centre
        rvAddWidthBase = rvWidth - vSeptumBaseRadialDisplacement + vSeptumThickness - lvFreeWallThickness + rvFreeWallThickness
        rvBaseJoinX = lvOuterRadius*math.cos(0.5*rvArcAroundBaseRadians)
        joinFactor = 2.0
        rvAddWidthRadiusBase = rvAddWidthBase - rvBaseJoinX - joinFactor*rvSideExtension

        # optimisation to convert RV inner height fraction to RV outer height using optimi
        lengthUpInner = 0.25*getApproximateEllipsePerimeter(lvInnerHeight, lvInnerRadius)
        rvInnerThetaUp = math.acos(rvInnerHeightFraction)
        rvThetaUp = 0.5*math.pi - (0.5*math.pi - rvInnerThetaUp)*elementsCountUpRV/(elementsCountUpRV - 0.5)
        for iter in range(100):
            lengthUpApex = getEllipseArcLength(lvInnerHeight, lvInnerRadius, 0.0, rvThetaUp)
            rvEstInnerThetaUp = updateEllipseAngleByArcLength(lvInnerHeight, lvInnerRadius, 0.0, lengthUpApex*(1.0 + 0.5/elementsCountUpLVApex))
            ratio = rvInnerThetaUp/rvEstInnerThetaUp
            rvThetaUp *= ratio
            if math.fabs(ratio - 1.0) < 1.0E-6:
                break
        else:
            print('Calculation of RV theta up did not converge after', iter + 1, 'iterations. Final ratio', ratio)
        rvOuterHeight = 0.5*lvOuterHeight*(1.0 + rvInnerHeightFraction)
        rvSideHeight = rvInnerHeightFraction*lvOuterHeight
        radialDisplacementStartRadiansUp = math.acos(rvOuterHeight/lvOuterHeight)

        # get LV inner points

        lvInnerx  = []
        lvInnerd1 = []
        lengthUpApex = getEllipseArcLength(lvInnerHeight, lvInnerRadius, 0.0, rvThetaUp)
        lengthUpRV = lengthUpInner - lengthUpApex
        elementSizeUpApex = lengthUpApex/elementsCountUpLVApex
        elementSizeUpRV = (lengthUpRV - 0.5*elementSizeUpApex)/(elementsCountUpRV - 0.5)
        elementSizeUpRVTransition = 0.5*(elementSizeUpApex + elementSizeUpRV)
        # LV apex points, noting s1, s2 is x, -y to get out outward s3
        lvApexInnerx = [ 0.0, 0.0, -lvInnerHeight ]
        lvApexInnerd1 = [ elementSizeUpApex, 0.0, 0.0 ]
        lvApexInnerd2 = [ 0.0, -elementSizeUpApex, 0.0, 0.0 ]
        lvApexInnerd3 = [ 0.0, 0.0, -lvApexThickness ]

        lvInnerRadiansUp = []
        radiansUp = 0.0
        for n2 in range(elementsCountUpLV):
            arcLengthUp = elementSizeUpApex if (n2 < elementsCountUpLVApex) else (elementSizeUpRVTransition if (n2 == elementsCountUpLVApex) else elementSizeUpRV)
            radiansUp = updateEllipseAngleByArcLength(lvInnerHeight, lvInnerRadius, radiansUp, arcLengthUp)
            if n2 == (elementsCountUpLV - 1):
                radiansUp = 0.5*math.pi
            lvInnerRadiansUp.append(radiansUp)

        for n2 in range(elementsCountUpLV):
            radiansUp = lvInnerRadiansUp[n2]
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            z = -lvInnerHeight*cosRadiansUp
            lvRadius = lvInnerRadius*sinRadiansUp
            rvArcAroundZRadians = rvArcAroundBaseRadians*(1.0 + (1.0 - rvArcApexFraction)*z/lvInnerHeight)

            # get radial displacement of centre of septum, a function of radiansUp
            xiUp = max(0.0, (radiansUp - radialDisplacementStartRadiansUp)/(0.5*math.pi - radialDisplacementStartRadiansUp))
            midSeptumDisplacement = interp.interpolateCubicHermite([0.0], [0.0], [vSeptumBaseRadialDisplacement], [0.0], xiUp)[0]

            nx, nd1 = getLeftVentricleInnerPoints(lvRadius, midSeptumDisplacement, rvArcAroundZRadians, z, \
                elementsCountAroundLVFreeWall, elementsCountAroundVSeptum, ivSulcusDerivativeFactor)

            lvInnerx.append(nx)
            lvInnerd1.append(nd1)

        # get ventricles outer points

        vOuterx  = []
        vOuterd1 = []
        lengthUpOuter = 0.25*getApproximateEllipsePerimeter(lvOuterHeight, lvOuterRadius)
        lengthUpApex = getEllipseArcLength(lvOuterHeight, lvOuterRadius, 0.0, rvThetaUp)
        lengthUpRV = lengthUpOuter - lengthUpApex
        elementSizeUpApex = lengthUpApex/elementsCountUpLVApex
        elementSizeUpRV = (lengthUpRV - 0.5*elementSizeUpApex)/(elementsCountUpRV - 0.5)
        elementSizeUpRVTransition = 0.5*(elementSizeUpApex + elementSizeUpRV)
        # apex points, noting s1, s2 is x, -y to get out outward s3
        lvApexOuterx = [ 0.0, 0.0, -lvOuterHeight ]
        lvApexOuterd1 = [ elementSizeUpApex, 0.0, 0.0 ]
        lvApexOuterd2 = [ 0.0, -elementSizeUpApex, 0.0, 0.0 ]
        lvApexOuterd3 = [ 0.0, 0.0, -lvApexThickness ]

        lvOuterRadiansUp = []
        radiansUp = 0.0
        for n2 in range(elementsCountUpLV):
            arcLengthUp = elementSizeUpApex if (n2 < elementsCountUpLVApex) else (elementSizeUpRVTransition if (n2 == elementsCountUpLVApex) else elementSizeUpRV)
            radiansUp = updateEllipseAngleByArcLength(lvOuterHeight, lvOuterRadius, radiansUp, arcLengthUp)
            if n2 == (elementsCountUpLV - 1):
                radiansUp = 0.5*math.pi
            lvOuterRadiansUp.append(radiansUp)

        for n2 in range(elementsCountUpLV):
            radiansUp = lvOuterRadiansUp[n2]
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            z = -lvOuterHeight*cosRadiansUp
            lvRadius = lvOuterRadius*sinRadiansUp

            rvArcAroundZRadians = rvArcAroundBaseRadians*(1.0 + (1.0 - rvArcApexFraction)*z/lvOuterHeight)

            # get size of RV at z
            xiUpWidth = 1.0 + z/rvOuterHeight
            xiUpSide = 1.0 + z/rvSideHeight
            widthExtension, sideExtension = getRVOuterSize(xiUpWidth, xiUpSide, rvAddWidthBase, rvWidthGrowthFactor, rvSideExtension, rvSideExtensionGrowthFactor)
            addWidthRadius = interp.interpolateHermiteLagrange([ 0.0 ], [ 0.0 ], [ rvAddWidthRadiusBase ], xiUpSide)[0]

            nx, nd1 = getVentriclesOuterPoints(lvRadius, widthExtension, sideExtension, addWidthRadius, rvArcAroundZRadians, z, \
                elementsCountAroundLVFreeWall, elementsCountAroundRVFreeWall, ivSulcusDerivativeFactor)

            vOuterx.append(nx)
            vOuterd1.append(nd1)

        # get radians around apex for each radial line on outer surface
        rvApexArcAroundRadians = rvArcAroundBaseRadians*rvArcApexFraction
        lvApexArcAroundRadians = 2.0*math.pi - rvApexArcAroundRadians
        lvApexRadiansPerElement = lvApexArcAroundRadians/(elementsCountAroundLVFreeWall + ivSulcusDerivativeFactor - 1.0)
        lvApexRadiansPerElementTransition = 0.5*lvApexRadiansPerElement*(1.0 + ivSulcusDerivativeFactor)
        ivSulcusRadiansPerElement = 2.0*lvApexRadiansPerElementTransition - lvApexRadiansPerElement
        rvApexRadiansPerElement = (rvApexArcAroundRadians - ivSulcusRadiansPerElement)/(elementsCountAroundRVFreeWall - 1)
        rvApexRadiansPerElementTransition = 0.5*(rvApexRadiansPerElement + ivSulcusRadiansPerElement)
        apexRadiansAround = []
        radiansAround = -0.5*rvApexArcAroundRadians
        for n1 in range(elementsCountAroundLVFreeWall):
            apexRadiansAround.append(radiansAround)
            if (n1 == 0) or (n1 == (elementsCountAroundLVFreeWall - 1)):
                radiansAround -= lvApexRadiansPerElementTransition
            else:
                radiansAround -= lvApexRadiansPerElement
        for n1 in range(elementsCountAroundRVFreeWall):
            apexRadiansAround.append(radiansAround)
            if (n1 == 0) or (n1 == (elementsCountAroundRVFreeWall - 1)):
                radiansAround -= rvApexRadiansPerElementTransition
            else:
                radiansAround -= rvApexRadiansPerElement

        # calculate derivative 2 on inner surfaces

        lvInnerd2 = []
        for n2 in range(elementsCountUpLV):
            lvInnerd2.append([])
        for n1 in range(pointsCountAroundOuter):
            cosRadiansAround = math.cos(apexRadiansAround[n1])
            sinRadiansAround = math.sin(apexRadiansAround[n1])
            apexd2 = [ (lvApexOuterd1[c]*cosRadiansAround + lvApexOuterd2[c]*sinRadiansAround) for c in range(3) ]
            nx = [ lvApexOuterx ]
            nd2 = [ apexd2 ]
            for n2 in range(elementsCountUpLV):
                nx.append(lvInnerx[n2][n1])
                nd2.append([ 0.0, 0.0, 0.0 ])
            nd2[-1] = [ 0.0, 0.0, elementSizeUpRV ]
            nd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True)
            for n2 in range(elementsCountUpLV):
                lvInnerd2[n2].append(nd2[n2 + 1])

        # calculate derivative 2 on outer surface

        vOuterd2 = []
        for n2 in range(elementsCountUpLV):
            vOuterd2.append([])
        for n1 in range(pointsCountAroundOuter):
            cosRadiansAround = math.cos(apexRadiansAround[n1])
            sinRadiansAround = math.sin(apexRadiansAround[n1])
            apexd2 = [ (lvApexOuterd1[c]*cosRadiansAround + lvApexOuterd2[c]*sinRadiansAround) for c in range(3) ]
            nx = [ lvApexOuterx ]
            nd2 = [ apexd2 ]
            for n2 in range(elementsCountUpLV):
                nx.append(vOuterx[n2][n1])
                nd2.append([ 0.0, 0.0, 0.0 ])
            nd2[-1] = [ 0.0, 0.0, elementSizeUpRV ]
            nd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True)
            for n2 in range(elementsCountUpLV):
                vOuterd2[n2].append(nd2[n2 + 1])

        # calculate derivative 3 on LV free wall and sub-RV apex from difference between inner and outer surfaces

        lvInnerd3 = []
        vOuterd3 = []
        for n2 in range(elementsCountUpLV):
            layerInnerd3 = []
            layerOuterd3 = []
            for n1 in range(pointsCountAroundOuter):
                if (n2 < elementsCountUpLVApex) or (n1 <= elementsCountAroundLVFreeWall):
                    ix = lvInnerx[n2][n1]
                    ox = vOuterx[n2][n1]
                    innerd3 = outerd3 = [ (ox[c] - ix[c]) for c in range(3) ]
                else:
                    innerd3 = outerd3 = [ 0.0, 0.0, 0.0 ]
                layerInnerd3.append(innerd3)
                layerOuterd3.append(outerd3)
            lvInnerd3.append(layerInnerd3)
            vOuterd3.append(layerOuterd3)

        # get points on inside of RV, anticlockwise starting at node opposite posterior interventricular sulcus

        septumOuterRadius = lvInnerRadius + vSeptumThickness
        rvInnerx = []
        rvInnerd1 = []
        rvInnerd2 = []
        rvInnerd3 = []
        for n2 in range(elementsCountUpLV):
            layerInnerx = []
            layerInnerd1 = []
            layerInnerd2 = []
            layerInnerd3 = []
            if n2 >= elementsCountUpLVApex:
                for n1 in range(elementsCountAroundLVFreeWall + 1, pointsCountAroundOuter):
                    outerx  = vOuterx [n2][n1]
                    outerd1 = vOuterd1[n2][n1]
                    outerd2 = vOuterd2[n2][n1]
                    unitNormal = vector.normalise(vector.crossproduct3(outerd1, outerd2))
                    innerd3 = vOuterd3[n2][n1] = vector.setMagnitude(unitNormal, rvFreeWallThickness)
                    innerx = [ (outerx[c] - innerd3[c]) for c in range(3) ]
                    # calculate inner d1 from curvature around
                    n1m = n1 - 1
                    n1p = (n1 + 1) % pointsCountAroundOuter
                    curvature = -0.5*(
                        interp.getCubicHermiteCurvature(vOuterx[n2][n1m], vOuterd1[n2][n1m], outerx, outerd1, unitNormal, 1.0) +
                        interp.getCubicHermiteCurvature(outerx, outerd1, vOuterx[n2][n1p], vOuterd1[n2][n1p], unitNormal, 0.0))
                    factor = 1.0 - curvature*rvFreeWallThickness
                    innerd1 = [ factor*c for c in outerd1 ]
                    # calculate inner d2 from curvature up
                    n2m = n2 - 1
                    curvature = -interp.getCubicHermiteCurvature(vOuterx[n2m][n1], vOuterd2[n2m][n1], outerx, outerd2, unitNormal, 1.0)
                    if n2 < (elementsCountUpLVApex - 1):
                        n2p = n2 + 1
                        curvature = 0.5*(curvature - interp.getCubicHermiteCurvature(outerx, outerd2, vOuterx[n2p][n1], vOuterd2[n2p][n1], unitNormal, 0.0))
                    factor = 1.0 - curvature*rvFreeWallThickness
                    innerd2 = [ factor*c for c in outerd2 ]
                    layerInnerx .append(innerx)
                    layerInnerd1.append(innerd1)
                    layerInnerd2.append(innerd2)
                    layerInnerd3.append(innerd3)
                # sample points onto outer septum
                innerRadiansUp = lvInnerRadiansUp[n2]
                z = -lvInnerHeight*math.cos(innerRadiansUp)
                radiansUp = math.acos(-z/lvOuterHeight)
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                lvRadius = septumOuterRadius*sinRadiansUp
                rvArcAroundZRadians = rvArcAroundBaseRadians*(1.0 + (1.0 - rvArcApexFraction)*z/lvOuterHeight)

                # get radial displacement of centre of septum, a function of radiansUp
                xiUp = max(0.0, (radiansUp - radialDisplacementStartRadiansUp)/(0.5*math.pi - radialDisplacementStartRadiansUp))
                midSeptumDisplacement = interp.interpolateCubicHermite([0.0], [0.0], [vSeptumBaseRadialDisplacement], [0.0], xiUp)[0]

                numberAround = 10
                nx, nd1 = getLeftVentricleInnerPoints(lvRadius, midSeptumDisplacement, rvArcAroundZRadians, z, numberAround, numberAround)
                # extract septum points, reversing order and derivative
                nx = [ layerInnerx [-1] ] + nx [-1:-numberAround:-1] + [ layerInnerx [1 - elementsCountAroundRVFreeWall] ]
                nd1 = [ layerInnerd1[-1] ] + [ ([-d for d in nd1[n1] ]) for n1 in range(-1, -numberAround, -1) ] + [ layerInnerd1[1 - elementsCountAroundRVFreeWall] ]
                scale = 2.0
                for n1 in [ -1, 0 ]:
                    nd1[n1] = [ scale*d for d in nd1[n1] ]
                px, pd1 = interp.sampleCubicHermiteCurves(nx, nd1, elementsCountAroundVSeptum + 2,
                    addLengthStart = (0.5/scale)*vector.magnitude(nd1[ 0]), lengthFractionStart = 0.5,
                    addLengthEnd = (0.5/scale)*vector.magnitude(nd1[-1]), lengthFractionEnd = 0.5,
                    arcLengthDerivatives = False)[0:2]
                for n1 in range(elementsCountAroundVSeptum + 1):
                    # d3 at ends of septum is toward adjacent interventricular sulcus, inside septum is across septum to lvInner
                    insideSeptum = (n1 > 0) and (n1 < elementsCountAroundVSeptum)
                    lvx = lvInnerx[n2][-n1] if insideSeptum else vOuterx[n2][-n1]
                    rvx = px[n1 + 1]
                    innerd3 = [ (lvx[c] - rvx[c]) for c in range(3) ]
                    layerInnerx .append(rvx)
                    layerInnerd1.append(pd1[n1 + 1])
                    layerInnerd2.append([ 0.0, 0.0, 0.0 ])
                    layerInnerd3.append(innerd3)
                    if insideSeptum:
                        lvInnerd3[n2][-n1] = [ -d for d in innerd3 ]
                # swizzle lists to start at node opposite posterior interventricular sulcus
                for li in [ layerInnerx, layerInnerd1, layerInnerd2, layerInnerd3 ]:
                    li.insert(0, li.pop())
            rvInnerx.append(layerInnerx)
            rvInnerd1.append(layerInnerd1)
            rvInnerd2.append(layerInnerd2)
            rvInnerd3.append(layerInnerd3)
        # calculate derivative 2 on inner RV septum
        n2Range = range(elementsCountUpLVApex, elementsCountUpLV)
        for n1 in range(-elementsCountAroundVSeptum, 1):
            nx = [ rvInnerx[n2][n1] for n2 in n2Range ]
            nd2 = [ rvInnerd2[n2][n1] for n2 in n2Range ]
            nd2 = interp.smoothCubicHermiteDerivativesLine([ rvInnerx[n2][n1] for n2 in n2Range ], [ rvInnerd2[n2][n1] for n2 in n2Range ])
            for n2 in n2Range:
                rvInnerd2[n2][n1] = nd2[n2 - elementsCountUpLVApex]

        # get points on RV apex curve, from posterior to anterior

        n2 = elementsCountUpLVApex
        dFactor = 2.0
        sx = []
        sd1 = []
        sd2 = []
        for n1 in range(1, elementsCountAroundRVFreeWall):
            r1 = 2*elementsCountAroundRVFreeWall - n1
            ax = rvInnerx[n2][r1]
            ad1 = [ -d for d in rvInnerd1[n2][r1] ]
            ad2 = [ -dFactor*d for d in rvInnerd2[n2][r1] ]
            bx = rvInnerx[n2][n1]
            bd1 = rvInnerd1[n2][n1]
            bd2 = [ dFactor*d for d in rvInnerd2[n2][n1] ]
            px, pd2, pe, pxi = interp.sampleCubicHermiteCurves([ ax, bx ], [ ad2, bd2 ], elementsCountOut = 2,
                addLengthStart = 0.5*vector.magnitude(ad2)/dFactor, lengthFractionStart = 0.5,
                addLengthEnd = 0.5*vector.magnitude(bd2)/dFactor, lengthFractionEnd = 0.5, arcLengthDerivatives = False)[0:4]
            sx .append(px[1])
            sd1.append(interp.interpolateSampleLinear([ ad1, bd1 ], pe[1:2], pxi[1:2])[0])
            sd2.append(pd2[1])
        ax  = rvInnerx [n2][0]
        ad1 = rvInnerd1[n2][0]
        ad2 = [ -d for d in rvInnerd2[n2][0] ]
        bx  = rvInnerx [n2][elementsCountAroundRVFreeWall]
        bd1 = [ -d for d in rvInnerd1[n2][elementsCountAroundRVFreeWall] ]
        bd2 = rvInnerd2[n2][elementsCountAroundRVFreeWall]
        px, pd1, pe, pxi = interp.sampleCubicHermiteCurves([ ax ] + sx + [ bx ], [ ad2 ] + sd1 + [ bd2 ], elementsCountAroundRVFreeWall + 2,
            addLengthStart = 0.5*vector.magnitude(ad2)/dFactor, lengthFractionStart = 0.5,
            addLengthEnd = 0.5*vector.magnitude(bd2)/dFactor, lengthFractionEnd = 0.5, arcLengthDerivatives = False)[0:4]
        pd2 = interp.interpolateSampleLinear([ ad1 ] + sd2 + [ bd1 ], pe, pxi)
        n2 = elementsCountUpLVApex - 1
        # loop skips first and last in sample:
        for n1 in range(1, elementsCountAroundRVFreeWall + 2):
            rvInnerx [n2].append(px [n1])
            rvInnerd1[n2].append(pd1[n1])
            o1 = (elementsCountAroundLVFreeWall + n1 - 1) % elementsCountAroundLV
            d1Factor = 1.0 if (n1 == 1) else (-1.0 if (n1 == (elementsCountAroundRVFreeWall + 1)) else 0.0)
            if (n1 == 1) or (n1 == elementsCountAroundRVFreeWall + 1):
                # collapsed RV corner uses triangle derivative d/dx{1|3} = +/-d2; outside d/dxi2 = +/-d1
                # compute derivative 2 to fit with glm vector on inside row:
                od2 = [ (d1Factor*lvInnerd1[n2][o1][c] + lvInnerd2[n2][o1][c] + lvInnerd3[n2][o1][c]) for c in range(3) ]
                id2 = interp.interpolateHermiteLagrangeDerivative(lvInnerx[n2][o1], od2, px[n1], 1.0)
                rvInnerd2[n2].append(id2)
            else:
                rvInnerd2[n2].append(pd2[n1])
            # compute derivative 3 to fit with glm vector on outside row:
            od3 = [ (d1Factor*vOuterd1[n2][o1][c] + vOuterd2[n2][o1][c] - vOuterd3[n2][o1][c]) for c in range(3) ]
            id3 = interp.interpolateHermiteLagrangeDerivative(vOuterx[n2][o1], od3, px[n1], 1.0)
            rvInnerd3[n2].append([ -d for d in id3 ])

        # create nodes on inner left ventricle

        apexNodeId = [ None, None ]

        node = nodes.createNode(nodeIdentifier, nodetemplateApex)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lvApexInnerx)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lvApexInnerd1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lvApexInnerd2)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lvApexInnerd3)
        apexNodeId[1] = nodeIdentifier
        nodeIdentifier = nodeIdentifier + 1

        lvInnerNodeId = []
        for n2 in range(elementsCountUpLV):
            nx  = lvInnerx[n2]
            nd1 = lvInnerd1[n2]
            nd2 = lvInnerd2[n2]
            nd3 = lvInnerd3[n2]
            layerNodeId = []
            for n1 in range(pointsCountAroundOuter):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nx[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, nd1[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, nd2[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, nd3[n1])
                layerNodeId.append(nodeIdentifier)
                nodeIdentifier = nodeIdentifier + 1
            lvInnerNodeId.append(layerNodeId)

        # create nodes on inner right ventricle

        rvInnerNodeId = []
        for n2 in range(elementsCountUpLV):
            nx  = rvInnerx[n2]
            nd1 = rvInnerd1[n2]
            nd2 = rvInnerd2[n2]
            nd3 = rvInnerd3[n2]
            layerNodeId = []
            for n1 in range(len(nx)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nx[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, nd1[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, nd2[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, nd3[n1])
                layerNodeId.append(nodeIdentifier)
                nodeIdentifier = nodeIdentifier + 1
            rvInnerNodeId.append(layerNodeId)
        # mirror RV apex so can index as for other rows
        n2 = elementsCountUpLVApex - 1
        for n1 in range(elementsCountAroundRVFreeWall - 1, 0, -1):
            rvInnerNodeId[n2].append(rvInnerNodeId[n2][n1])

        # create nodes on outer ventricles

        apexNodeId = [ None, None ]

        node = nodes.createNode(nodeIdentifier, nodetemplateApex)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lvApexOuterx)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lvApexOuterd1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lvApexOuterd2)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lvApexOuterd3)
        apexNodeId[1] = nodeIdentifier
        nodeIdentifier = nodeIdentifier + 1

        vOuterNodeId = []
        for n2 in range(elementsCountUpLV):
            nx  = vOuterx[n2]
            nd1 = vOuterd1[n2]
            nd2 = vOuterd2[n2]
            nd3 = vOuterd3[n2]
            layerNodeId = []
            for n1 in range(pointsCountAroundOuter):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nx[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, nd1[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, nd2[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, nd3[n1])
                layerNodeId.append(nodeIdentifier)
                nodeIdentifier = nodeIdentifier + 1
            vOuterNodeId.append(layerNodeId)

        #################
        # Create elements
        #################

        lvMeshGroup = lvGroup.getMeshGroup(mesh)
        rvMeshGroup = rvGroup.getMeshGroup(mesh)
        vSeptumMeshGroup = vSeptumGroup.getMeshGroup(mesh)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = tricubichermite.createEftNoCrossDerivatives()

        elementIdentifier = 1

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        norl = elementsCountAroundLV
        nowl = 2 + elementsCountAroundRVFreeWall + elementsCountAroundLV*elementsCountUpLV + elementsCountAroundVSeptum*2*elementsCountUpRV
        norr = 2*elementsCountAroundRVFreeWall

        for e2 in range(elementsCountUpLV):

            if e2 < elementsCountUpLVApex:

                # LV apex

                meshGroups = [ lvMeshGroup ]
                for e1 in range(elementsCountAroundLV):

                    eft1 = eft
                    scalefactors = None

                    va = e1
                    vb = (va + 1)%elementsCountAroundLV
                    if e2 == 0:
                        # create apex elements, varying eft scale factor identifiers around apex
                        # scale factor identifiers follow convention of offsetting by 100 for each 'version'
                        nids = [ 1       , 2 + va       , 2 + vb,
                                 1 + nowl, 2 + va + nowl, 2 + vb + nowl ]
                        eft1 = tricubichermite.createEftShellPoleBottom(va*100, vb*100)
                        # calculate general linear map coefficients
                        radiansa = apexRadiansAround[va]
                        radiansb = apexRadiansAround[vb]
                        dRadiansa = 0.5*(apexRadiansAround[va - 1] - apexRadiansAround[vb])
                        if dRadiansa < 0.0:
                            dRadiansa += math.pi
                        dRadiansb = 0.5*(apexRadiansAround[va] - apexRadiansAround[va + 2 - elementsCountAroundLV])
                        if dRadiansb < 0.0:
                            dRadiansb += math.pi
                        scalefactors = [ -1.0,
                            math.cos(radiansa), math.sin(radiansa), dRadiansa,
                            math.cos(radiansb), math.sin(radiansb), dRadiansb,
                            math.cos(radiansa), math.sin(radiansa), dRadiansa,
                            math.cos(radiansb), math.sin(radiansb), dRadiansb
                        ]
                    else:
                        bni1 = 2 + norl*(e2 - 1)
                        nids = [ bni1        + va, bni1        + vb, bni1        + norl + va, bni1        + norl + vb,
                                 bni1 + nowl + va, bni1 + nowl + vb, bni1 + nowl + norl + va, bni1 + nowl + norl + vb ]

                    result1 = elementtemplate1.defineField(coordinates, -1, eft1)
                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    result3 = element.setScaleFactors(eft1, scalefactors) if scalefactors else None
                    #print('create element lv apex', elementIdentifier, result1, result2, result3, nids)
                    elementIdentifier = elementIdentifier + 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

            if e2 >= elementsCountUpLVApex:

                # LV free wall elements starting from -1 = anterior interventricular sulcus collapsed element

                for e1 in range(-1, elementsCountAroundLVFreeWall):

                    eft1 = eft
                    scalefactors = None
                    meshGroups = [ lvMeshGroup ]

                    va = e1 if (e1 >= 0) else (elementsCountAroundLV - 1)
                    vb = e1 + 1
                    nids = [ lvInnerNodeId[e2 - 1][va], lvInnerNodeId[e2 - 1][vb], lvInnerNodeId[e2][va], lvInnerNodeId[e2][vb],
                             vOuterNodeId [e2 - 1][va], vOuterNodeId [e2 - 1][vb], vOuterNodeId [e2][va], vOuterNodeId [e2][vb] ]
                    if e1 == -1:
                        # anterior interventricular sulcus, collapsed to 6 node wedge
                        nids[0] = rvInnerNodeId[e2 - 1][elementsCountAroundRVFreeWall]
                        nids[2] = rvInnerNodeId[e2    ][elementsCountAroundRVFreeWall]
                        nids.pop(6)
                        nids.pop(4)
                        meshGroups += [ rvMeshGroup ]
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                        if e2 == elementsCountUpLVApex:
                            # collapsed RV corner uses triangle derivative d/dx1 = -d2; outside d/dxi2 = d1
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        ln_map = [ 1, 2, 3, 4, 5, 5, 6, 6 ]
                        remapEftLocalNodes(eft1, 6, ln_map)

                    result1 = elementtemplate1.defineField(coordinates, -1, eft1)

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    result3 = element.setScaleFactors(eft1, scalefactors) if scalefactors else None
                    #print('create element lv', elementIdentifier, result1, result2, result3, nids)
                    elementIdentifier = elementIdentifier + 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

            if e2 >= (elementsCountUpLVApex - 1):

                # RV free wall elements starting from -1 = posterior interventricular sulcus collapsed element

                for e1 in range(-1, elementsCountAroundRVFreeWall):

                    eft1 = eft
                    scalefactors = None
                    meshGroups = [ rvMeshGroup ]
                    ua = e1
                    ub = e1 + 1
                    va = elementsCountAroundLVFreeWall + e1
                    vb = (va + 1)%elementsCountAroundLV
                    e2m = max(e2 - 1, elementsCountUpLVApex - 1)
                    nids = [ rvInnerNodeId[e2m][ua], rvInnerNodeId[e2m][ub], rvInnerNodeId[e2][ua], rvInnerNodeId[e2][ub],
                             vOuterNodeId [e2m][va], vOuterNodeId [e2m][vb], vOuterNodeId [e2][va], vOuterNodeId [e2][vb] ]
                    if e2 == (elementsCountUpLVApex - 1):
                        if e1 == -1:
                            continue
                        nids[0] = lvInnerNodeId[e2m][va]
                        nids[1] = lvInnerNodeId[e2m][vb]
                        nids.pop(7)
                        nids.pop(6)
                        meshGroups += [ lvMeshGroup ]
                        # collapsed elements at RV apex
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS2, [])
                        if e1 == 0:
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                            remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        elif e1 == (elementsCountAroundVSeptum - 1):
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                            remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 3, 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                            remapEftNodeValueLabel(eft1, [ 7, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        ln_map = [ 1, 2, 3, 4, 5, 6, 5, 6 ]
                        remapEftLocalNodes(eft1, 6, ln_map)
                    elif e1 == -1:
                        # posterior interventricular sulcus, collapsed to 6 node wedge
                        nids[0] = lvInnerNodeId[e2 - 1][elementsCountAroundLVFreeWall]
                        nids[2] = lvInnerNodeId[e2    ][elementsCountAroundLVFreeWall]
                        nids.pop(6)
                        nids.pop(4)
                        meshGroups += [ lvMeshGroup ]
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                        if e2 == elementsCountUpLVApex:
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            # collapsed RV corner uses triangle derivative d/dx1 = d2; outside d/dxi2 = -d1
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        ln_map = [ 1, 2, 3, 4, 5, 5, 6, 6 ]
                        remapEftLocalNodes(eft1, 6, ln_map)
                    elif e1 == 0:
                        # general linear map d3 adjacent to collapsed posterior interventricular sulcus
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        if e2 == elementsCountUpLVApex:
                            # collapsed RV corner uses outside d/dxi2 = -d1
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                    elif e1 == (elementsCountAroundRVFreeWall - 1):
                        # general linear map d3 adjacent to collapsed anterior interventricular sulcus
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        if e2 == elementsCountUpLVApex:
                            # collapsed RV corner uses outside d/dxi2 = d1
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                    elif e2 == elementsCountUpLVApex:
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])

                    result1 = elementtemplate1.defineField(coordinates, -1, eft1)

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    result3 = element.setScaleFactors(eft1, scalefactors) if scalefactors else None
                    #print('create element rv', elementIdentifier, result1, result2, result3, nids)
                    elementIdentifier = elementIdentifier + 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

            # interventricular septum

            if e2 >= elementsCountUpLVApex:

                for e1 in range(elementsCountAroundVSeptum):

                    eft1 = eft
                    scalefactors = None
                    meshGroups = [ lvMeshGroup, rvMeshGroup, vSeptumMeshGroup ]

                    ua = -e1
                    ub = ua - 1
                    va = elementsCountAroundLVFreeWall + e1
                    vb = (va + 1)%elementsCountAroundLV
                    nids = [ lvInnerNodeId[e2 - 1][va], lvInnerNodeId[e2 - 1][vb], lvInnerNodeId[e2][va], lvInnerNodeId[e2][vb],
                             rvInnerNodeId[e2 - 1][ua], rvInnerNodeId[e2 - 1][ub], rvInnerNodeId[e2][ua], rvInnerNodeId[e2][ub] ]
                    if e2 == elementsCountUpLVApex:
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        if e1 == 0:
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            # general linear map d3 adjacent to collapsed posterior interventricular sulcus
                            remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            # collapsed RV corner uses triangle derivative d/dx3 = d2; outside d/dxi2 = -d1
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            scaleEftNodeValueLabels(eft1, [ 7 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                            remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            scaleEftNodeValueLabels(eft1, [ 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                        elif e1 == (elementsCountAroundVSeptum - 1):
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            # collapsed RV corner uses triangle derivative d/dx3 = d2; outside d/dxi2 = d1
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                            # general linear map d3 adjacent to collapsed anterior interventricular sulcus
                            remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            scaleEftNodeValueLabels(eft1, [ 7 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                            scaleEftNodeValueLabels(eft1, [ 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                            remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            scaleEftNodeValueLabels(eft1, [ 7, 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                    else:
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        if e1 == 0:
                            # general linear map d3 adjacent to collapsed posterior interventricular sulcus
                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                            remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            scaleEftNodeValueLabels(eft1, [ 6, 8 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                        elif e1 == (elementsCountAroundRVFreeWall - 1):
                            # general linear map d3 adjacent to collapsed anterior interventricular sulcus
                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                            scaleEftNodeValueLabels(eft1, [ 5, 7 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                            remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        else:
                            scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])

                    result1 = elementtemplate1.defineField(coordinates, -1, eft1)

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    result3 = element.setScaleFactors(eft1, scalefactors) if scalefactors else None
                    #print('create element septum', elementIdentifier, result1, result2, result3, nids)
                    elementIdentifier = elementIdentifier + 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

        # apex annotation point
        element1 = mesh.findElementByIdentifier(1)
        markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
        nodeIdentifier += 1
        cache.setNode(markerPoint)
        markerName.assignString(cache, 'APEX')  # interlex.org has 'apex of heart'
        markerLocation.assignMeshLocation(cache, element1, [ 0.0, 0.0, 1.0 ])

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        Stops at end of ventricles, hence can be called from ventriclesbase.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        elementsCountAroundLVFreeWall = options['Number of elements around LV free wall']
        elementsCountAroundRVFreeWall = options['Number of elements around RV free wall']
        elementsCountAroundLV = elementsCountAroundLVFreeWall + elementsCountAroundRVFreeWall
        elementsCountAroundVSeptum = elementsCountAroundRVFreeWall
        elementsCountUpLVApex = options['Number of elements up LV apex']
        elementsCountUpRV = options['Number of elements up RV']
        elementsCountUpLV = elementsCountUpLVApex + elementsCountUpRV
        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughLVWall = options['Refine number of elements through LV wall']
        refineElementsCountThroughRVWall = options['Refine number of elements through wall']

        startPostApexElementIdentifier = elementsCountAroundLV*elementsCountUpLVApex + 1
        elementsCountPostApexLayer = elementsCountAroundLV + 2 + elementsCountAroundVSeptum
        lastVentriclesElementIdentifier = startPostApexElementIdentifier + elementsCountAroundVSeptum + elementsCountUpRV*elementsCountPostApexLayer - 1

        element = meshrefinement._sourceElementiterator.next()
        while element.isValid():
            numberInXi1 = refineElementsCountSurface
            numberInXi2 = refineElementsCountSurface
            numberInXi3 = refineElementsCountThroughLVWall
            elementIdentifier = element.getIdentifier()
            if elementIdentifier >= startPostApexElementIdentifier:
                n1 = (elementIdentifier - startPostApexElementIdentifier - elementsCountAroundVSeptum)
                if n1 < 0:
                    # collapsed elements at RV apex
                    numberInXi2 = refineElementsCountThroughLVWall
                else:
                    n1 = n1 % elementsCountPostApexLayer
                    if (n1 == 0) or (n1 == (elementsCountAroundLVFreeWall + 1)):
                        # collapsed elements on posterior or anterior interventricular sulcus:
                        numberInXi1 = refineElementsCountThroughLVWall
                        #print(n1,'refine collapse element', elementIdentifier, numberInXi1, numberInXi2, numberInXi3)
                    elif (n1 > (elementsCountAroundLVFreeWall + 1)) and (n1 <= (elementsCountAroundLV + 1)):
                        numberInXi3 = refineElementsCountThroughRVWall
            meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
            if elementIdentifier == lastVentriclesElementIdentifier:
                return  # finish on last so can continue in ventriclesbase
            element = meshrefinement._sourceElementiterator.next()


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
        # create endocardium and epicardium groups
        fm = region.getFieldmodule()
        lvGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("left ventricle myocardium"))
        rvGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("right ventricle myocardium"))
        vSeptumGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("interventricular septum"))
        mesh2d = fm.findMeshByDimension(2)
        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi3_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))
        is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_lv = lvGroup.getFieldElementGroup(mesh2d)
        is_rv = rvGroup.getFieldElementGroup(mesh2d)
        is_lv_endo = fm.createFieldAnd(is_lv, is_exterior_face_xi3_0)
        is_rv_endo = fm.createFieldOr(fm.createFieldAnd(fm.createFieldAnd(is_rv, is_exterior_face_xi3_0),
                                                        fm.createFieldNot(is_lv_endo)),
                                      fm.createFieldAnd(vSeptumGroup.getFieldElementGroup(mesh2d), is_exterior_face_xi3_1))
        is_v_epi = fm.createFieldAnd(fm.createFieldOr(is_lv, is_rv),
                                     fm.createFieldAnd(is_exterior_face_xi3_1,
                                                     fm.createFieldNot(vSeptumGroup.getFieldElementGroup(mesh2d))))
        epiGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_heart_term("epicardium"))
        epiGroup.getMeshGroup(mesh2d).addElementsConditional(is_v_epi)
        lvEndoGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_heart_term("endocardium of left ventricle"))
        lvEndoGroup.getMeshGroup(mesh2d).addElementsConditional(is_lv_endo)
        rvEndoGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_heart_term("endocardium of right ventricle"))
        rvEndoGroup.getMeshGroup(mesh2d).addElementsConditional(is_rv_endo)


def getSeptumPoints(septumArcRadians, lvRadius, radialDisplacement, elementsCountAroundLVFreeWall, elementsCountAroundVSeptum, z, n3):
    '''
    Symmetric 2-cubic interpolation of n septum elements around arc.
    :return: sx[], sd1[]
    '''
    radiansPerElementAroundLVFreeWall = (2.0*math.pi - septumArcRadians)/elementsCountAroundLVFreeWall
    # get cubic curve with arc length scaling across to centre of septum
    radiansAround = 0.5*septumArcRadians
    circleArcLength = lvRadius*radiansAround
    x1 = [ lvRadius - radialDisplacement, 0.0, z ]
    d1 = [ 0.0, circleArcLength, 0.0 ]
    x2 = [ lvRadius*math.cos(radiansAround), lvRadius*math.sin(radiansAround), z ]
    d2 = [ -circleArcLength*math.sin(radiansAround), circleArcLength*math.cos(radiansAround), 0.0 ]
    cubicLength = interp.computeCubicHermiteArcLength(x1, d1, x2, d2, False)
    elementLengthMid = (2.0*cubicLength - lvRadius*radiansPerElementAroundLVFreeWall)/(elementsCountAroundVSeptum + n3 - 1)
    odd = elementsCountAroundVSeptum % 2
    sx, sd1 = interp.sampleCubicHermiteCurves([ x1, x2 ], [ d1, d2 ], elementsCountAroundVSeptum//2 + odd,
        lengthFractionStart = 0.5 if odd else 1.0,
        addLengthEnd = cubicLength - 0.5*elementsCountAroundVSeptum*elementLengthMid, arcLengthDerivatives = True)[0:2]
    if odd:
        sx = sx[1:]
        sd1 = sd1[1:]
    return sx, sd1


def getLeftVentricleInnerPoints(lvRadius, midSeptumDisplacement, septumArcAroundRadians, z,
    elementsCountAroundLVFreeWall, elementsCountAroundVSeptum, ivSulcusDerivativeFactor = 1.0):
    '''
    Get array of points and derivatives around inside of leftventricle anticlockwise starting from
    anterior interventricular junction.
    LV is assumed to be circular, centred at x, y = (0,0), but optionally displaced inward over septum.
    :param lvRadius: Radius of the LV.
    :param midSeptumDisplacement: Inward displacement of mid septum.
    :param septumArcAroundRadians: Angle in radians around septum.
    :param z: Coordinate to give to all x[2]. All d1[2] are zero.
    :param elementsCountAroundLVFreeWall: Number of elements around LV part of outer wall.
    :param elementsCountAroundVSeptum: Number of elements around Ventricular septum.
    :param ivSulcusDerivativeFactor: Factor, typically less than 1.0 giving derivative magnitude
    at interventricular sulcus as a fraction of LV derivative, to allow a sharper bend.
    :return: Arrays nx[], nd[].
    '''
    lvArcAroundRadians = 2.0*math.pi - septumArcAroundRadians
    lvLength = lvArcAroundRadians*lvRadius
    elementSizeLV = lvLength/(elementsCountAroundLVFreeWall + ivSulcusDerivativeFactor - 1.0)
    elementSizeLVTransition = 0.5*elementSizeLV*(1.0 + ivSulcusDerivativeFactor)
    a = lvRadius - midSeptumDisplacement
    b = lvRadius
    circleLimit = 0.5*math.pi*lvRadius
    ellipseLength = 0.25*getApproximateEllipsePerimeter(a, b)
    sepLength = 2.0*(circleLimit + ellipseLength) - lvLength
    ivSulcusDerivative = ivSulcusDerivativeFactor*elementSizeLV
    elementSizeSep = (sepLength - ivSulcusDerivative)/(elementsCountAroundVSeptum - 1)
    elementSizeSepTransition = 0.5*(elementSizeSep + ivSulcusDerivative)
    # get 2-D points half way around, from mid LV free wall to mid septum
    hx = []
    hd = []
    lvOdd = (elementsCountAroundLVFreeWall % 2)
    sepOdd = (elementsCountAroundVSeptum % 2)
    nIvSulcus = elementsCountAroundLVFreeWall//2
    distance = 0.5*elementSizeLV if lvOdd else 0.0
    halfPointsCount = nIvSulcus + 1 + elementsCountAroundVSeptum//2
    for n in range(halfPointsCount):
        # get coordinate x and derivative direction d (magnitude is set later)
        if distance < circleLimit:
            angleRadians = -math.pi + distance/lvRadius
            x = [ lvRadius*math.cos(angleRadians), lvRadius*math.sin(angleRadians) ]
            d = [ -math.sin(angleRadians), math.cos(angleRadians) ]
        else:
            angleRadians = updateEllipseAngleByArcLength(a, b, -0.5*math.pi, distance - circleLimit)
            x = [ a*math.cos(angleRadians), b*math.sin(angleRadians) ]
            d = [ -a*math.sin(angleRadians), b*math.cos(angleRadians) ]
        if n < nIvSulcus:
            magnitude = elementSizeLV
            distance += elementSizeLV if (n < (nIvSulcus - 1)) else elementSizeLVTransition
        elif n == nIvSulcus:
            magnitude = ivSulcusDerivative
            distance += elementSizeSepTransition
        else:
            magnitude = elementSizeSep
            distance += elementSizeSep
        hx.append(x)
        hd.append(vector.setMagnitude(d, magnitude))
    #distance -= (0.5*elementSizeSep if sepOdd else elementSizeSep)

    nx = []
    nd = []
    lvMirrorPointsLimit = nIvSulcus + lvOdd
    for p in range(lvMirrorPointsLimit):
        n = nIvSulcus - p
        x = hx[n]
        nx.append([ x[0], -x[1], z ])
        d = hd[n]
        nd.append([ -d[0], d[1], 0.0 ])
    for n in range(halfPointsCount):
        x = hx[n]
        nx.append([ x[0], x[1], z ])
        d = hd[n]
        nd.append([ d[0], d[1], 0.0 ])
    sepMirrorPointsStart = -1 if sepOdd else -2
    sepMirrorPointsLimit = -(elementsCountAroundVSeptum//2 + 1)
    for p in range(sepMirrorPointsStart, sepMirrorPointsLimit, -1):
        n = halfPointsCount + p
        x = hx[n]
        nx.append([ x[0], -x[1], z ])
        d = hd[n]
        nd.append([ -d[0], d[1], 0.0 ])
    return nx, nd


def getVentriclesOuterPoints(lvRadius, widthExtension, sideExtension, addWidthRadius, rvArcAroundRadians, z,
        elementsCountAroundLVFreeWall, elementsCountAroundRVFreeWall, ivSulcusDerivativeFactor = 1.0):
    '''
    Get array of points and derivatives around outside of ventricles anticlockwise starting from
    anterior interventricular sulcus.
    LV is assumed to be circular, centred at x, y = (0,0)
    :param lvRadius: Radius of the LV.
    :param widthExtension: Amount to add to LV to get RV width at centre.
    :param sideExtension: Additional radius to add only laterally around septum.
    :param addWidthRadius: Amount to add to width radius without affecting RV width, varying
        from zero with zero derivative to max at base.
    :param rvArcAroundRadians: Angle in radians around outside of RV to tangent nodes where it meets LV.
    :param z: Coordinate to give to all x[2]. All d1[2] are zero.
    :param elementsCountAroundLVFreeWall: Number of elements around LV part of outer wall.
    :param elementsCountAroundRVFreeWall: Number of elements around RV part of outer wall.
    :param ivSulcusDerivativeFactor: Factor, typically less than 1.0 giving derivative magnitude
    at interventricular sulcus as a fraction of LV derivative, to allow a sharper bend.
    :return: Arrays nx[], nd[].
    '''
    lvArcAroundRadians = 2.0*math.pi - rvArcAroundRadians
    lvLength = lvArcAroundRadians*lvRadius
    elementSizeLV = lvLength/(elementsCountAroundLVFreeWall + ivSulcusDerivativeFactor - 1.0)
    elementSizeLVTransition = 0.5*elementSizeLV*(1.0 + ivSulcusDerivativeFactor)
    a = lvRadius + addWidthRadius  # width radius
    b = lvRadius + sideExtension   # side radius
    ellipseCentrex = lvRadius + widthExtension - a
    joinRadians = -0.5*max(math.pi, rvArcAroundRadians)
    x1 = [ lvRadius*math.cos(joinRadians), lvRadius*math.sin(joinRadians) ]
    d1 = [ -math.sin(joinRadians), math.cos(joinRadians) ]
    x2 = [ ellipseCentrex, -b ]
    d2 = [ 1.0, 0.0 ]
    cubicLength = interp.computeCubicHermiteArcLength(x1, d1, x2, d2, True)
    d1 = [ d1[0]*cubicLength, d1[1]*cubicLength ]
    d2 = [ d2[0]*cubicLength, d2[1]*cubicLength ]
    circleLimit = lvRadius*(math.pi + joinRadians)
    cubicLimit = circleLimit + cubicLength
    ellipseLength = 0.25*getApproximateEllipsePerimeter(a, b)
    rvLength = 2.0*(circleLimit + cubicLength + ellipseLength) - lvLength
    ivSulcusDerivative = ivSulcusDerivativeFactor*elementSizeLV
    elementSizeRV = (rvLength - ivSulcusDerivative)/(elementsCountAroundRVFreeWall - 1)
    elementSizeRVTransition = 0.5*(elementSizeRV + ivSulcusDerivative)
    # get 2-D points half way around, from mid LV free wall to mid RV free wall
    hx = []
    hd = []
    lvOdd = (elementsCountAroundLVFreeWall % 2)
    rvOdd = (elementsCountAroundRVFreeWall % 2)
    nIvSulcus = elementsCountAroundLVFreeWall//2
    distance = 0.5*elementSizeLV if lvOdd else 0.0
    halfPointsCount = nIvSulcus + 1 + elementsCountAroundRVFreeWall//2
    for n in range(halfPointsCount):
        # get coordinate x and derivative direction d (magnitude is set later)
        if distance < circleLimit:
            angleRadians = -math.pi + distance/lvRadius
            x = [ lvRadius*math.cos(angleRadians), lvRadius*math.sin(angleRadians) ]
            d = [ -math.sin(angleRadians), math.cos(angleRadians) ]
        elif distance < cubicLimit:
            xi = (distance - circleLimit)/cubicLength
            x = interp.interpolateCubicHermite(x1, d1, x2, d2, xi)
            d = interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, xi)
        else:
            angleRadians = updateEllipseAngleByArcLength(a, b, -0.5*math.pi, distance - cubicLimit)
            x = [ ellipseCentrex + a*math.cos(angleRadians), b*math.sin(angleRadians) ]
            d = [ -a*math.sin(angleRadians), b*math.cos(angleRadians) ]
        if n < nIvSulcus:
            magnitude = elementSizeLV
            distance += elementSizeLV if (n < (nIvSulcus - 1)) else elementSizeLVTransition
        elif n == nIvSulcus:
            magnitude = ivSulcusDerivative
            distance += elementSizeRVTransition
        else:
            magnitude = elementSizeRV
            distance += elementSizeRV
        hx.append(x)
        hd.append(vector.setMagnitude(d, magnitude))
    #distance -= (0.5*elementSizeRV if rvOdd else elementSizeRV)

    nx = []
    nd = []
    lvMirrorPointsLimit = nIvSulcus + lvOdd
    for p in range(lvMirrorPointsLimit):
        n = nIvSulcus - p
        x = hx[n]
        nx.append([ x[0], -x[1], z ])
        d = hd[n]
        nd.append([ -d[0], d[1], 0.0 ])
    for n in range(halfPointsCount):
        x = hx[n]
        nx.append([ x[0], x[1], z ])
        d = hd[n]
        nd.append([ d[0], d[1], 0.0 ])
    rvMirrorPointsStart = -1 if rvOdd else -2
    rvMirrorPointsLimit = -(elementsCountAroundRVFreeWall//2 + 1)
    for p in range(rvMirrorPointsStart, rvMirrorPointsLimit, -1):
        n = halfPointsCount + p
        x = hx[n]
        nx.append([ x[0], -x[1], z ])
        d = hd[n]
        nd.append([ -d[0], d[1], 0.0 ])
    return nx, nd


def getRVOuterSize(xiUpWidth, xiUpSide, rvWidth, rvWidthGrowthFactor, rvSideExtension, rvSideExtensionGrowthFactor):
    '''
    :return: widthExtension, sideExtension
    '''
    if xiUpWidth < 0.0:
        return 0.0, 0.0
    xiUpWidthRaw = interp.getCubicHermiteCurvesPointAtArcDistance([ [ 0.0 ], [ 1.0 ] ],
        [ [ 2.0*(1.0 - rvWidthGrowthFactor) ] , [ 2.0*rvWidthGrowthFactor ] ],
        xiUpWidth)[3]
    xiUpWidthMod = 1.0 - (1.0 - xiUpWidthRaw)*(1.0 - xiUpWidthRaw)
    xiUpSideMod = interp.getCubicHermiteCurvesPointAtArcDistance([ [ 0.0 ], [ 1.0 ] ],
        [ [ 2.0*(1.0 - rvSideExtensionGrowthFactor) ] , [ 2.0*rvSideExtensionGrowthFactor ] ],
        max(xiUpSide, 0.0))[3]
    widthExtension = interp.interpolateCubicHermite([0.0], [0.0], [rvWidth], [0.0], xiUpWidthMod)[0]
    sideExtension = interp.interpolateCubicHermite([0.0], [0.0], [rvSideExtension], [0.0], xiUpSideMod)[0]
    return widthExtension, sideExtension
