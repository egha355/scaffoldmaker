"""
Generates 3-D Left and Right ventricles mesh starting from modified sphere shell mesh.
"""

from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.annotation.heart_terms import get_heart_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, setEftScaleFactorIds
from scaffoldmaker.utils.geometry import getApproximateEllipsePerimeter, getEllipseArcLength, updateEllipseAngleByArcLength
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class MeshType_3d_heartventricles2(Scaffold_base):
    '''
    classdocs
    '''

    @staticmethod
    def getName():
        return '3D Heart Ventricles 2'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements around LV free wall' : 5,
            'Number of elements around ventricular septum' : 7,
            'Number of elements up LV apex' : 1,
            'Number of elements up ventricular septum' : 4,
            'LV outer height' : 0.95,
            'LV outer radius' : 0.5,
            'LV free wall thickness' : 0.12,
            'LV apex thickness' : 0.06,
            'RV inner height' : 0.8,
            'RV arc around degrees' : 200.0,
            'RV free wall thickness' : 0.04,
            'RV width' : 0.4,
            'RV extra cross radius base' : 0.1,
            'Ventricular septum thickness' : 0.1,
            'Ventricular septum base radial displacement' : 0.15,
            'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements surface' : 4,
            'Refine number of elements through LV wall' : 1,
            'Refine number of elements through RV wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around LV free wall',
            'Number of elements around ventricular septum',
            'Number of elements up LV apex',
            'Number of elements up ventricular septum',
            'LV outer height',
            'LV outer radius',
            'LV free wall thickness',
            'LV apex thickness',
            'RV inner height',
            'RV arc around degrees',
            'RV free wall thickness',
            'RV width',
            'RV extra cross radius base',
            'Ventricular septum thickness',
            'Ventricular septum base radial displacement',
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through LV wall',
            'Refine number of elements through RV wall'
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
            'Refine number of elements through RV wall',
            'Number of elements up LV apex']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Number of elements around LV free wall',
            'Number of elements up ventricular septum']:
            if options[key] < 2:
                options[key] = 2
        for key in [
            'Number of elements around ventricular septum']:
            if options[key] < 3:
                options[key] = 3
        for key in [
            'LV outer height',
            'LV outer radius',
            'LV free wall thickness',
            'LV apex thickness',
            'RV inner height',
            'RV free wall thickness',
            'RV width',
            'RV extra cross radius base',
            'Ventricular septum thickness',
            'Ventricular septum base radial displacement']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['RV inner height'] > 0.99*options['LV outer height']:
            options['RV inner height'] = 0.99*options['LV outer height']
        if options['RV arc around degrees'] < 45.0:
            options['RV arc around degrees'] = 45.0
        elif options['RV arc around degrees'] > 270.0:
            options['RV arc around degrees'] = 270.0
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
        elementsCountAroundVSeptum = options['Number of elements around ventricular septum']
        elementsCountAroundLV = elementsCountAroundLVFreeWall + elementsCountAroundVSeptum
        elementsCountUpLVApex = options['Number of elements up LV apex']
        elementsCountUpVSeptum = options['Number of elements up ventricular septum']
        elementsCountUpLV = elementsCountUpLVApex + elementsCountUpVSeptum
        elementsCountUpRV = elementsCountUpVSeptum + 1
        elementsCountAroundRV = elementsCountAroundVSeptum + 2
        lvOuterHeight = options['LV outer height']
        lvOuterRadius = options['LV outer radius']
        lvFreeWallThickness = options['LV free wall thickness']
        lvApexThickness = options['LV apex thickness']
        rvInnerHeight = options['RV inner height']
        rvArcAroundRadians = math.radians(options['RV arc around degrees'])
        rvFreeWallThickness = options['RV free wall thickness']
        rvWidth = options['RV width']
        rvExtraCrossRadiusBase = options['RV extra cross radius base']
        vSeptumThickness = options['Ventricular septum thickness']
        vSeptumBaseRadialDisplacement = options['Ventricular septum base radial displacement']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)
        cache = fm.createFieldcache()

        lvGroup = AnnotationGroup(region, get_heart_term("left ventricle myocardium"))
        rvGroup = AnnotationGroup(region, get_heart_term("right ventricle myocardium"))
        vSeptumGroup = AnnotationGroup(region, get_heart_term("interventricular septum"))
        annotationGroups = [ lvGroup, rvGroup, vSeptumGroup ]

        #################
        # Create nodes
        #################

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        nodetemplateApex = nodetemplate

        nodeIdentifier = 1

        # LV nodes

        rvOuterWidthBase = rvWidth - vSeptumBaseRadialDisplacement + vSeptumThickness - lvFreeWallThickness + rvFreeWallThickness

        # get element spacing up
        # approximate correction to move up to next inner node, since bottom of RV is halfway up element
        rvHeightCorrection = (elementsCountUpVSeptum - 1)/(elementsCountUpVSeptum - 0.5)
        rvThetaUp = math.acos(rvInnerHeight*rvHeightCorrection/lvOuterHeight)
        rvOuterHeight = 0.5*(rvInnerHeight + lvOuterHeight)
        rvCrossHeight = rvInnerHeight
        radialDisplacementStartRadiansUp = math.acos(rvOuterHeight/lvOuterHeight)
        lvArcAroundRadians = 2.0*math.pi - rvArcAroundRadians

        norl = elementsCountAroundLV
        nowl = elementsCountUpLV*elementsCountAroundLV + 1
        for n3 in range(2):
            a = lvOuterHeight
            b = lvOuterRadius
            bVS = lvOuterRadius - lvFreeWallThickness + vSeptumThickness
            if n3 == 0:
                a -= lvApexThickness
                b -= lvFreeWallThickness
                bVS = b

            lengthUp = 0.25*getApproximateEllipsePerimeter(a, b)
            lengthUpApex = getEllipseArcLength(a, b, 0.0, rvThetaUp)
            lengthUpSeptum = lengthUp - lengthUpApex
            elementSizeUpSeptum = lengthUpSeptum/(elementsCountUpVSeptum - 1)
            #if n3 == 1:
            #    elementSizeUpApex = (lengthUpApex - (2.0/3.0)*elementSizeUpSeptum)/(elementsCountUpLVApex + (2.0/3.0))
            #    elementSizeUpSeptumTransition = (2.0/3.0)*(elementSizeUpApex + elementSizeUpSeptum)
            #else:
            elementSizeUpApex = (lengthUpApex - 0.5*elementSizeUpSeptum)/(elementsCountUpLVApex + 0.5)
            elementSizeUpSeptumTransition = 0.5*(elementSizeUpApex + elementSizeUpSeptum)
            # have reduced derivative around RV junction to fit extra row of elements:
            d2FactorRvTransition = 0.5*elementSizeUpSeptumTransition/elementSizeUpApex  #  2.0/3.0

            #print('length total apex septum', lengthUp, lengthUpApex, lengthUpSeptum, ' size apex septum', elementSizeUpApex, elementSizeUpSeptum)

            # apex node, noting s1, s2 is x, -y to get out outward s3
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ 0.0, 0.0, -a ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ elementSizeUpApex, 0.0, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ 0.0, -elementSizeUpApex, 0.0, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, -lvApexThickness ])
            nodeIdentifier = nodeIdentifier + 1

            radiansUp = 0.0
            for n2 in range(elementsCountUpLV):

                if n2 < elementsCountUpLVApex:
                    arcLengthUp = elementSizeUpApex
                    dArcLengthUp = elementSizeUpApex
                else:
                    if n2 == elementsCountUpLVApex:
                        arcLengthUp = elementSizeUpSeptumTransition
                    else:
                        arcLengthUp = elementSizeUpSeptum
                    dArcLengthUp = elementSizeUpSeptum
                radiansUp = updateEllipseAngleByArcLength(a, b, radiansUp, arcLengthUp)
                if n2 == (elementsCountUpLV - 1):
                    radiansUp = 0.5*math.pi
                #print(n3, n2, ' -> ', math.degrees(radiansUp), 'dArcLengthUp', dArcLengthUp)
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                z = -a*cosRadiansUp
                radius = b*sinRadiansUp
                bSeptum = b
                if (n3 == 1) and (n2 >= elementsCountUpLVApex):
                    bSeptum -= (lvFreeWallThickness - vSeptumThickness)
                radiusSeptum = bSeptum*sinRadiansUp

                elementSizeCrossSeptum = 0.5*(radius*lvArcAroundRadians/elementsCountAroundLVFreeWall + vSeptumThickness*sinRadiansUp)
                elementSizeAroundLVFreeWall = (radius*lvArcAroundRadians - elementSizeCrossSeptum)/(elementsCountAroundLVFreeWall - 1)
                elementSizeAroundLVFreeWallTransition = 0.5*(elementSizeAroundLVFreeWall + elementSizeCrossSeptum)

                # get radial displacement of centre of septum, a function of radiansUp
                xiUp = max(0.0, (radiansUp - radialDisplacementStartRadiansUp)/(0.5*math.pi - radialDisplacementStartRadiansUp))
                radialDisplacement = interp.interpolateCubicHermite([0.0], [0.0], [vSeptumBaseRadialDisplacement], [0.0], xiUp)[0]

                radiansAround = -0.5*rvArcAroundRadians
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                if (n3 == 1) and (n2 < elementsCountUpLVApex):
                    # Follow RV on outside below septum
                    dEndMag = elementSizeCrossSeptum
                    xiUpZ = 1.0 + z/rvOuterHeight
                    xiUpCross = 1.0 + z/rvCrossHeight
                    rvAddWidthRadius, rvAddCrossRadius = getRVOuterSize(xiUpZ, xiUpCross, rvOuterWidthBase, rvExtraCrossRadiusBase)
                    sx, sd1 = getRvOuterPoints(rvArcAroundRadians, radius, rvAddWidthRadius, rvAddCrossRadius, elementsCountAroundVSeptum, dEndMag, z, xiUpCross)
                else:
                    sx, sd1 = getSeptumPoints(rvArcAroundRadians, radiusSeptum, radialDisplacement, elementsCountAroundLVFreeWall, elementsCountAroundVSeptum, z, n3)

                # do finite difference to get d/dxi2
                dxi = 1.0E-3
                ds_dxi = dArcLengthUp
                if (n3 == 1) and (n2 < elementsCountUpLVApex):
                    dzr_dRadiansUp = [ a*sinRadiansUp, b*cosRadiansUp ]
                else:
                    dzr_dRadiansUp = [ a*sinRadiansUp, bSeptum*cosRadiansUp ]
                ds_dRadiansUp = vector.magnitude(dzr_dRadiansUp)
                dzr_ds = vector.normalise(dzr_dRadiansUp)
                dz = dxi*ds_dxi*dzr_ds[0]
                dr = dxi*ds_dxi*dzr_ds[1]
                if (n3 == 1) and (n2 < elementsCountUpLVApex):
                    # Follow RV on outside below septum
                    radiansUpPlus = radiansUp + math.sqrt(dz*dz + dr*dr)/ds_dRadiansUp
                    dEndMag = 0.5*((radius + dr)*lvArcAroundRadians/elementsCountAroundLVFreeWall + vSeptumThickness*math.sin(radiansUpPlus))
                    xiUpZPlus = 1.0 + (z + dz)/rvOuterHeight
                    xiUpCross = 1.0 + (z + dz)/rvCrossHeight
                    rvAddWidthRadius, rvAddCrossRadius = getRVOuterSize(xiUpZPlus, xiUpCross, rvOuterWidthBase, rvExtraCrossRadiusBase)
                    tx, td1 = getRvOuterPoints(rvArcAroundRadians, radius + dr, rvAddWidthRadius, rvAddCrossRadius, elementsCountAroundVSeptum, dEndMag, z + dz, xiUpCross)
                else:
                    dxiUp_dxi = ds_dxi/(ds_dRadiansUp*(0.5*math.pi - radialDisplacementStartRadiansUp))
                    dRadialDisplacement = dxi*dxiUp_dxi*interp.interpolateCubicHermiteDerivative([0.0], [0.0], [vSeptumBaseRadialDisplacement], [0.0], xiUp)[0]
                    tx, td1 = getSeptumPoints(rvArcAroundRadians, radiusSeptum + dr, radialDisplacement + dRadialDisplacement, elementsCountAroundLVFreeWall, elementsCountAroundVSeptum, z + dz, n3)
                # true values for LV freewall
                dzr_dRadiansUp = [ a*sinRadiansUp, b*cosRadiansUp ]
                dzr_ds = vector.normalise(dzr_dRadiansUp)

                for n1 in range(elementsCountAroundLV):

                    if (n1 > 0) and (n1 < elementsCountAroundVSeptum):
                        n1mid = (elementsCountAroundVSeptum + 1)//2
                        if n1 < n1mid:
                            nr = elementsCountAroundVSeptum//2 - n1
                            x = [ sx[nr][0], -sx[nr][1], sx[nr][2] ]
                            dx_ds1 = [ -sd1[nr][0], sd1[nr][1], sd1[nr][2] ]
                            xt = [ tx[nr][0], -tx[nr][1], tx[nr][2] ]
                        else:
                            nr = n1 - n1mid
                            x = sx[nr]
                            dx_ds1 = sd1[nr]
                            xt = tx[nr]
                        dx_ds2 = [ (xt[c] - x[c])/dxi for c in range(3) ]
                    else:
                        if n1 == 0:
                            radiansAround = -0.5*rvArcAroundRadians
                            dMagAround = elementSizeCrossSeptum
                        elif n1 == elementsCountAroundVSeptum:
                            radiansAround = 0.5*rvArcAroundRadians
                            dMagAround = elementSizeCrossSeptum
                        else:
                            radiansAround = 0.5*rvArcAroundRadians + (elementSizeAroundLVFreeWallTransition + (n1 - elementsCountAroundVSeptum - 1)*elementSizeAroundLVFreeWall)/radius
                            dMagAround = elementSizeAroundLVFreeWall
                        cosRadiansAround = math.cos(radiansAround)
                        sinRadiansAround = math.sin(radiansAround)
                        x = [ radius*cosRadiansAround, radius*sinRadiansAround, z ]
                        dx_ds1 = [ dMagAround*-sinRadiansAround, dMagAround*cosRadiansAround, 0.0 ]
                        dx_ds2 = [ cosRadiansAround*ds_dxi*dzr_ds[1], sinRadiansAround*ds_dxi*dzr_ds[1], ds_dxi*dzr_ds[0] ]

                    if (n3 == 1) and (n2 == (elementsCountUpLVApex - 1)) and (n1 >= 0) and (n1 <= elementsCountAroundVSeptum):
                        dx_ds2 = [ d*d2FactorRvTransition for d in dx_ds2 ]

                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    if n3 == 1:
                        # get derivative3 from difference in coordinates
                        innerNode = nodes.findNodeByIdentifier(nodeIdentifier - nowl)
                        cache.setNode(innerNode)
                        result, xInner = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                        dx_ds3 = [ (x[c] - xInner[c]) for c in range(3) ]
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    nodeIdentifier = nodeIdentifier + 1

        # RV nodes

        nidr = nodeIdentifier

        rxOuter = []
        rd1Outer = []
        rd2Outer = []
        rd3Outer = []
        elementSizeUpRvTransition1 = d2FactorRvTransition*elementSizeUpApex
        elementSizeUpRv = (lengthUp - elementsCountUpLVApex*elementSizeUpApex - 1.5*elementSizeUpRvTransition1)/(elementsCountUpRV - 1.5)
        elementSizeUpRvTransition2 = 0.5*(elementSizeUpRvTransition1 + elementSizeUpRv)
        distUp = elementSizeUpApex*elementsCountUpLVApex + elementSizeUpRvTransition1
        radiansUp = updateEllipseAngleByArcLength(a, b, 0.0, distUp)

        #print('RV element size apex', elementSizeUpApex, ', transition1', elementSizeUpRvTransition2, ', ', transition2', elementSizeUpRvTransition1, ', rv', elementSizeUpRv, \
        #    ', total', elementsCountUpLVApex*elementsCountUpLVApex + elementSizeUpRvTransition1 + elementSizeUpRvTransition2 + (elementsCountUpRV - 2)*elementSizeUpRv, ' vs ', lengthUp)

        for n2 in range(elementsCountUpRV):
            if n2 == (elementsCountUpRV - 1):
                radiansUp = 0.5*math.pi
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)

            radiansAround = -0.5*rvArcAroundRadians
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)

            lvOuterRadius = b*math.sin(radiansUp)
            ax = [ lvOuterRadius*cosRadiansAround, lvOuterRadius*sinRadiansAround, -a*cosRadiansUp ]
            z = ax[2]
            elementSizeCrossSeptum = 0.5*(lvOuterRadius*lvArcAroundRadians/elementsCountAroundLVFreeWall + vSeptumThickness*sinRadiansUp)
            dEndMag = elementSizeCrossSeptum
            ad2 = [ b*cosRadiansUp*cosRadiansAround, b*cosRadiansUp*sinRadiansAround, a*sinRadiansUp ]
            scale = elementSizeUpRv/vector.magnitude(ad2)
            if n2 == 0:
                scale *= elementSizeUpRvTransition1/elementSizeUpRv
            ad2 = [ d*scale for d in ad2 ]

            xiUpZ =1.0 + z/rvOuterHeight
            xiUpCross = 1.0 + z/rvCrossHeight

            rvAddWidthRadius, rvAddCrossRadius = getRVOuterSize(xiUpZ, xiUpCross, rvOuterWidthBase, rvExtraCrossRadiusBase)
            #print('z',z,'xiUpZ', xiUpZ, 'xiUpFast', xiUpFast, 'xiUpSlow', xiUpSlow)
            #print('rvAddWidthRadius =', rvAddWidthRadius, ', rvAddCrossRadius=', rvAddCrossRadius)
            rx, rd1 = getRvOuterPoints(rvArcAroundRadians, lvOuterRadius, rvAddWidthRadius, rvAddCrossRadius, elementsCountAroundRV, dEndMag, z, xiUpCross, True)

            if n2 < (elementsCountUpRV - 1):
                # get points at node position plus increment of xi2, to calculate derivative 2
                dxi = 1.0E-3
                xPlus = [ ax[i] + ad2[i]*dxi for i in range(3) ]
                lvOuterRadiusPlus = math.sqrt(xPlus[0]*xPlus[0] + xPlus[1]*xPlus[1])
                dArcLengthUp = elementSizeUpRv
                ds_dxi = dArcLengthUp
                dzr_dRadiansUp = [ a*sinRadiansUp, b*cosRadiansUp ]
                ds_dRadiansUp = vector.magnitude(dzr_dRadiansUp)
                dzr_ds = vector.normalise(dzr_dRadiansUp)
                dz = dxi*ds_dxi*dzr_ds[0]
                dr = dxi*ds_dxi*dzr_ds[1]
                radiansUpPlus = radiansUp + math.sqrt(dz*dz + dr*dr)/ds_dRadiansUp
                dEndMagPlus = 0.5*(lvOuterRadiusPlus*lvArcAroundRadians/elementsCountAroundLVFreeWall + vSeptumThickness*math.sin(radiansUpPlus))
                xiUpZPlus = 1.0 + xPlus[2]/rvOuterHeight
                xiUpCross = 1.0 + xPlus[2]/rvCrossHeight
                rvAddWidthRadiusPlus, rvAddCrossRadiusPlus = getRVOuterSize(xiUpZPlus, xiUpCross, rvOuterWidthBase, rvExtraCrossRadiusBase)
                px, pd1 = getRvOuterPoints(rvArcAroundRadians, lvOuterRadiusPlus, rvAddWidthRadiusPlus, rvAddCrossRadiusPlus, elementsCountAroundRV, dEndMagPlus, xPlus[2], xiUpCross, True)

                rd2 = [ [(px[n][c] - rx[n][c])/dxi for c in range(3)] for n in range(len(rx)) ]
            else:
                rd2 = [ad2]*len(rx)

            #if n2 == 0:
            #    for n in range(len(rd2)):
            #        rd2[n] = [ d*elementSizeUpRvTransition1/elementSizeUpRv for d in rd2[n] ]

            rd3 = []
            for n in range(len(rx)):
                normal = vector.crossproduct3(rd1[n], rd2[n])
                mag = rvFreeWallThickness/vector.magnitude(normal)
                rd3.append([ c*mag for c in normal ])

            rxOuter.append(rx)
            rd1Outer.append(rd1)
            rd2Outer.append(rd2)
            rd3Outer.append(rd3)
            radiansUp = updateEllipseAngleByArcLength(a, b, radiansUp, elementSizeUpRv if (n2 > 0) else elementSizeUpRvTransition2)

        # get inner RV nodes from outer

        nidl = 2 + (elementsCountUpLVApex - 1)*norl

        rxInner = []
        rd1Inner = []
        rd2Inner = []
        rd3Inner = []
        # fix bottom row below
        for n2 in range(1, elementsCountUpRV):
            ix = []
            id1 = []
            id2 = []
            id3 = []
            rx = rxOuter[n2]
            rd1 = rd1Outer[n2]
            rd2 = rd2Outer[n2]
            rd3 = rd3Outer[n2]

            for n1 in range(1, elementsCountAroundRV - 2):
                ix.append([ (rx[n1][c] - rd3[n1][c]) for c in range(3) ] )
                unitRadial = vector.normalise(rd3[n1])

                # calculate inner d1 from curvature around
                curvature = 0.0
                count = 0
                if n1 > 0:
                    curvature -= interp.getCubicHermiteCurvature(rx[n1 - 1], rd1[n1 - 1], rx[n1], rd1[n1], unitRadial, 1.0)
                    count += 1
                if n1 < (elementsCountAroundRV - 2):
                    curvature -= interp.getCubicHermiteCurvature(rx[n1], rd1[n1], rx[n1 + 1], rd1[n1 + 1], unitRadial, 0.0)
                    count += 1
                curvature /= count
                factor = 1.0 - curvature*rvFreeWallThickness
                id1.append([ factor*c for c in rd1[n1] ])

                # calculate inner d2 from curvature up
                curvature = 0.0
                count = 0
                if n2 > 0:
                    curvature -= interp.getCubicHermiteCurvature(rxOuter[n2 - 1][n1], rd2Outer[n2 - 1][n1], rx[n1], rd2[n1], unitRadial, 1.0)
                    count += 1
                if n2 < (elementsCountUpRV - 1):
                    curvature -= interp.getCubicHermiteCurvature(rx[n1], rd2[n1], rxOuter[n2 + 1][n1], rd2Outer[n2 + 1][n1], unitRadial, 0.0)
                    count += 1
                curvature /= count
                factor = 1.0 - curvature*rvFreeWallThickness
                id2.append([ factor*c for c in rd2Outer[n2][n1] ])

                id3.append(rd3[n1])

            # RV inlet/outlet radius

            for n1 in [ 0, elementsCountAroundRV - 2 ]:
                # interpolate node inside RV side edge
                n1r = 0 if (n1 == 0) else (elementsCountAroundRV - 3)
                nid = nidl + nowl + norl*n2 + (1 if (n1 == 0) else (elementsCountAroundVSeptum - 1))
                node = nodes.findNodeByIdentifier(nid)
                cache.setNode(node)
                result, ax = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                result, ad1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                result, ad2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                ad1 = [ -d for d in ad1 ]
                bx = ix[n1r]
                bd1 = id1[n1r]
                bd2 = id2[n1r]
                ad1 = [ 2.0*d for d in ad1 ]
                bd1 = [ 2.0*d for d in bd1 ]
                if n1 > 0:
                    ax, ad1, ad2, bx, bd1, bd2 = bx, bd1, bd2, ax, ad1, ad2
                x = interp.interpolateCubicHermite(ax, ad1, bx, bd1, 0.5)
                dx_ds1 = interp.interpolateCubicHermiteDerivative(ax, ad1, bx, bd1, 0.5)
                dx_ds1 = [ 0.5*d for d in dx_ds1 ]
                dx_ds2 = [ 0.5*(ad2[c] + bd2[c]) for c in range(3) ]
                ox = rxOuter[n2][n1]
                dx_ds3 = [ (ox[c] - x[c]) for c in range(3) ]
                ix.insert(n1, x)
                id1.insert(n1, dx_ds1)
                id2.insert(n1, dx_ds2)
                id3.insert(n1, dx_ds3)
                rd3Outer[n2][n1] = dx_ds3

            rxInner.append(ix)
            rd1Inner.append(id1)
            rd2Inner.append(id2)
            rd3Inner.append(id3)

        # Bottom of RV
        ix = []
        id1 = []
        id2 = []
        id3 = []
        for n1 in range(1, elementsCountAroundRV - 2):
            # interpolate node inside RV bottom edge
            nid = nidl + nowl + norl + n1
            #print(n2, 'n1', n1, 'n1r', n1r, 'LV nid', nid, )
            node = nodes.findNodeByIdentifier(nid)
            cache.setNode(node)
            result, ax = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, ad1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            result, ad2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
            ad2 = [ -d for d in ad2 ]
            bx = rxInner[0][n1]
            bd1 = rd1Inner[0][n1]
            bd2 = rd2Inner[0][n1]
            ad2 = [ 1.5*d for d in ad2 ]
            bd2 = [ 1.5*d for d in bd2 ]
            x = interp.interpolateCubicHermite(ax, ad2, bx, bd2, 0.5)
            dx_ds1 = [ 0.5*(ad1[c] + bd1[c]) for c in range(3) ]
            dx_ds2 = interp.interpolateCubicHermiteDerivative(ax, ad2, bx, bd2, 0.5)
            dx_ds2 = [ (1.0/3.0)*d for d in dx_ds2 ]
            ox = rxOuter[0][n1]
            dx_ds3 = [ (ox[c] - x[c]) for c in range(3) ]
            ix.append(x)
            id1.append(dx_ds1)
            id2.append(dx_ds2)
            id3.append(dx_ds3)
            rd3Outer[0][n1] = dx_ds3

        # Bottom corners of RV
        for n1 in [ 0, elementsCountAroundRV - 2 ]:
            # interpolate node inside RV side edge
            n1r = 0 if (n1 == 0) else (elementsCountAroundRV - 2)
            ax = rxInner[0][n1r]
            ad1 = rd2Inner[0][n1r]
            n1r = 0 if (n1 == 0) else (elementsCountAroundRV - 3)
            bx = ix[n1r]
            bd1 = id1[n1r]
            if n1 == 0:
                ad1 = [ -d for d in ad1 ]
            else:
                ax, ad1, bx, bd1 = bx, bd1, ax, ad1
            ad1 = [ 1.5*d for d in ad1 ]
            bd1 = [ 1.5*d for d in bd1 ]
            x = interp.interpolateCubicHermite(ax, ad1, bx, bd1, 0.5)
            dx_ds1 = interp.interpolateCubicHermiteDerivative(ax, ad1, bx, bd1, 0.5)
            dx_ds1 = [ (1.0/3.0)*d for d in dx_ds1 ]
            if n1 == 0:
                dx_ds2 = [ -d for d in dx_ds1 ]
            else:
                dx_ds2 = dx_ds1
            ox = rxOuter[0][n1]
            dx_ds3 = [ (ox[c] - x[c]) for c in range(3) ]
            ix.insert(n1, x)
            id1.insert(n1, dx_ds1)
            id2.insert(n1, dx_ds2)
            id3.insert(n1, dx_ds3)
            rd3Outer[0][n1] = dx_ds3

        rxInner.insert(0, ix)
        rd1Inner.insert(0, id1)
        rd2Inner.insert(0, id2)
        rd3Inner.insert(0, id3)

        zero = [ 0.0, 0.0, 0.0 ]
        rv_nids = [[],[]]
        for n3 in range(2):
            for n2 in range(elementsCountUpRV):
                nids = []
                rv_nids[n3].append(nids)
                if n3 == 0:
                    rx = rxInner[n2]
                    rd1 = rd1Inner[n2]
                    rd2 = rd2Inner[n2]
                    rd3 = rd3Inner[n2]
                else:
                    rx = rxOuter[n2]
                    rd1 = rd1Outer[n2]
                    rd2 = rd2Outer[n2]
                    rd3 = rd3Outer[n2]
                for n1 in range(elementsCountAroundRV - 1):
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rx[n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rd1[n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rd2[n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, rd3[n1])
                    nid = nodeIdentifier
                    nodeIdentifier += 1
                    nids.append(nid)


        #################
        # Create elements
        #################

        mesh = fm.findMeshByDimension(3)

        lvMeshGroup = lvGroup.getMeshGroup(mesh)
        rvMeshGroup = rvGroup.getMeshGroup(mesh)
        vSeptumMeshGroup = vSeptumGroup.getMeshGroup(mesh)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = tricubichermite.createEftNoCrossDerivatives()

        elementIdentifier = 1

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        norr = elementsCountAroundRV - 1
        nowr = norr*elementsCountUpRV

        # LV elements

        # calculate values used for apex scale factors
        radiansPerElementAroundLVFreeWall = lvArcAroundRadians/elementsCountAroundLVFreeWall
        if elementsCountAroundVSeptum > 2:
            radiansPerElementAroundSeptum = (rvArcAroundRadians - radiansPerElementAroundLVFreeWall)/(elementsCountAroundVSeptum - 1)
            radiansPerElementAroundSeptumTransition = 0.5*(radiansPerElementAroundLVFreeWall + radiansPerElementAroundSeptum)
        else:
            radiansPerElementAroundSeptum = radiansPerElementAroundSeptumTransition = rvArcAroundRadians/elementsCountAroundVSeptum
        radiansAround = -0.5*rvArcAroundRadians

        for e2 in range(elementsCountUpLV):

            for e1 in range(elementsCountAroundLV):

                eft1 = eft
                scalefactors = None
                meshGroups = [ lvMeshGroup ]

                if e2 == 0:
                    # create bottom apex elements, varying eft scale factor identifiers around apex
                    # scale factor identifiers follow convention of offsetting by 100 for each 'version'
                    va = e1
                    vb = (e1 + 1)%elementsCountAroundLV
                    nids = [ 1       , 2 + va       , 2 + vb,
                             1 + nowl, 2 + va + nowl, 2 + vb + nowl ]
                    eft1 = tricubichermite.createEftShellPoleBottom(va*100, vb*100)
                    # calculate general linear map coefficients
                    if e1 < elementsCountAroundVSeptum:
                        deltaRadians = dRadians1 = dRadians2 = radiansPerElementAroundSeptum
                        if e1 == 0:
                            deltaRadians = radiansPerElementAroundSeptumTransition
                            dRadians1 = radiansPerElementAroundLVFreeWall
                        elif e1 == (elementsCountAroundVSeptum - 1):
                            deltaRadians = radiansPerElementAroundSeptumTransition
                            dRadians2 = radiansPerElementAroundLVFreeWall
                    else:
                        deltaRadians = dRadians1 = dradians2 = radiansPerElementAroundLVFreeWall
                    radiansAroundNext = radiansAround + deltaRadians
                    radians1 = -radiansAround
                    radians2 = -radiansAroundNext
                    scalefactors = [ -1.0,
                        math.cos(radians1), math.sin(radians1), dRadians1,
                        math.cos(radians2), math.sin(radians2), dRadians2,
                        math.cos(radians1), math.sin(radians1), dRadians1,
                        math.cos(radians2), math.sin(radians2), dRadians2
                    ]
                    radiansAround = radiansAroundNext
                else:
                    bnil = 2 + norl*(e2 - 1) + e1
                    bnjl = 2 + norl*(e2 - 1) + (e1 + 1)%elementsCountAroundLV
                    nids = [ bnil       , bnjl       , bnil + norl       , bnjl + norl, \
                             bnil + nowl, bnjl + nowl, bnil + norl + nowl, bnjl + norl + nowl ]
                    if (e2 >= elementsCountUpLVApex) and (e1 < elementsCountAroundVSeptum):
                        meshGroups += [ rvMeshGroup, vSeptumMeshGroup ]
                        if (e2 == elementsCountUpLVApex) or (e1 == 0):
                            nids[4] = nidr + norr*(e2 - elementsCountUpLVApex) + e1
                        if e1 == 0:
                            nids[6] = nids[4] + norr
                        if (e2 == elementsCountUpLVApex) or (e1 == (elementsCountAroundVSeptum - 1)):
                            nids[5] = nidr + norr*(e2 - elementsCountUpLVApex) + e1 + 1
                        if e1 == (elementsCountAroundVSeptum - 1):
                            nids[7] = nids[5] + norr
                        if e2 == elementsCountUpLVApex:
                            if e1 == 0:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                scalefactors = [ -1.0 ]
                            elif e1 == (elementsCountAroundVSeptum - 1):
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                scalefactors = [ -1.0 ]
                            else:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                scalefactors = [ -1.0 ]
                        else:
                            if e1 == 0:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                scalefactors = [ -1.0 ]
                            elif e1 == (elementsCountAroundVSeptum - 1):
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                scalefactors = [ -1.0 ]

                result = elementtemplate1.defineField(coordinates, -1, eft1)

                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                result3 = element.setScaleFactors(eft1, scalefactors) if scalefactors else None
                #print('create element lv', elementIdentifier, result, result2, result3, nids)
                elementIdentifier = elementIdentifier + 1

                for meshGroup in meshGroups:
                    meshGroup.addElement(element)

        # RV elements

        nidl = 2 + (elementsCountUpLVApex - 1)*norl

        for e2 in range(elementsCountUpRV):

            for e1 in range(elementsCountAroundRV):

                eft1 = eft
                scalefactors = None
                meshGroups = [ rvMeshGroup ]

                if e2 == 0:
                    meshGroups += [ lvMeshGroup ]
                    if (e1 == 0) or (e1 == (elementsCountAroundRV - 1)):
                        # skip as two fewer elements on bottom row
                        continue
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    scalefactors = [ -1.0 ]
                    if e1 == 1:
                        remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    elif e1 == (elementsCountAroundRV - 2):
                        remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    else:
                        remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 3, 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                    bnil = nidl + norl*e2 + e1 - 1
                    bnjl = bnil + 1
                    bnir = nidr + norr*e2 + e1 - 1
                    bnjr = bnir + 1
                    nids = [ bnil       , bnjl       , bnir       , bnjr, \
                             bnil + nowl, bnjl + nowl, bnir + nowr, bnjr + nowr ]
                else:
                    if e1 == 0:
                        meshGroups += [ lvMeshGroup ]
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        if e2 == 1:
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        bnil = nidl + norl*(e2 - 1) + e1
                        bnir = nidr + norr*(e2 - 1) + e1
                        nids = [ bnil       , bnir       , bnil + norl       , bnir + norr, \
                                 bnil + nowl, bnir + nowr, bnil + norl + nowl, bnir + norr + nowr ]
                    elif e1 == (elementsCountAroundRV - 1):
                        meshGroups += [ lvMeshGroup ]
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        if e2 == 1:
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        bnil = nidl + norl*(e2 - 1) + elementsCountAroundVSeptum
                        bnir = nidr + norr*(e2 - 1) + elementsCountAroundRV - 2
                        nids = [ bnir       , bnil       , bnir + norr       , bnil + norl, \
                                 bnir + nowr, bnil + nowl, bnir + norr + nowr, bnil + norl + nowl ]
                        scalefactors = [ -1.0 ]
                    else:
                        bnir = nidr + norr*(e2 - 1) + e1 - 1
                        bnjr = bnir + 1
                        nids = [ bnir       , bnjr       , bnir + norr       , bnjr + norr, \
                                 bnir + nowr, bnjr + nowr, bnir + norr + nowr, bnjr + norr + nowr ]

                result = elementtemplate1.defineField(coordinates, -1, eft1)

                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                result3 = element.setScaleFactors(eft1, scalefactors) if scalefactors else None
                #print('create element rv', elementIdentifier, result, result2, result3, nids)
                elementIdentifier = elementIdentifier + 1

                for meshGroup in meshGroups:
                    meshGroup.addElement(element)

        fm.endChange()
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
        elementsCountAroundVSeptum = options['Number of elements around ventricular septum']
        elementsCountAroundLV = elementsCountAroundLVFreeWall + elementsCountAroundVSeptum
        elementsCountUpLVApex = options['Number of elements up LV apex']
        elementsCountUpVSeptum = options['Number of elements up ventricular septum']
        elementsCountUpLV = elementsCountUpLVApex + elementsCountUpVSeptum
        elementsCountUpRV = elementsCountUpVSeptum + 1
        elementsCountAroundRV = elementsCountAroundVSeptum + 2
        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughLVWall = options['Refine number of elements through LV wall']
        refineElementsCountThroughRVWall = options['Refine number of elements through RV wall']

        startRvElementIdentifier = elementsCountAroundLV*elementsCountUpLV + 1
        limitRvElementIdentifier = startRvElementIdentifier + elementsCountUpRV*elementsCountAroundRV - 2

        element = meshrefinement._sourceElementiterator.next()
        while element.isValid():
            numberInXi1 = refineElementsCountSurface
            numberInXi2 = refineElementsCountSurface
            elementId = element.getIdentifier()
            if elementId < startRvElementIdentifier:
                numberInXi3 = refineElementsCountThroughLVWall
            elif elementId < limitRvElementIdentifier:
                numberInXi3 = refineElementsCountThroughRVWall
            meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
            if elementId == (limitRvElementIdentifier - 1):
                return  # finish on last so can continue in ventriclesbase
            element = meshrefinement._sourceElementiterator.next()

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup for mesh.
        """
        if not options['Refine']:
            return cls.generateBaseMesh(region, options)
        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)
        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        cls.refineMesh(meshrefinement, options)
        return meshrefinement.getAnnotationGroups()


def getSeptumPoints(septumArcRadians, lvRadius, radialDisplacement, elementsCountAroundLVFreeWall, elementsCountAroundVSeptum, z, n3):
        '''
        Symmetric 2-cubic interpolation of n septum elements around arc.
        :return: x[], dx_ds1[]
        '''
        radiansPerElementAroundLVFreeWall = (2.0*math.pi - septumArcRadians)/elementsCountAroundLVFreeWall
        # get cubic curve with arc length scaling across to centre of septum
        radiansAround = 0.5*septumArcRadians
        circleArcLength = lvRadius*radiansAround
        v1 = [ lvRadius - radialDisplacement, 0.0, z ]
        d1 = [ 0.0, circleArcLength, 0.0 ]
        v2 = [ lvRadius*math.cos(radiansAround), lvRadius*math.sin(radiansAround), z ]
        d2 = [ -circleArcLength*math.sin(radiansAround), circleArcLength*math.cos(radiansAround), 0.0 ]
        cubicArcLength = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, False)
        scale = cubicArcLength/circleArcLength
        d1 = [ d*scale for d in d1 ]
        d2 = [ d*scale for d in d2 ]
        elementLengthMid = (2.0*cubicArcLength - lvRadius*radiansPerElementAroundLVFreeWall)/(elementsCountAroundVSeptum + n3 - 1)
        #elementLengthMid = (2.0*cubicArcLength - lvRadius*radiansPerElementAroundLVFreeWall)/(elementsCountAroundVSeptum + n3*0.5 - 1)
        #elementLengthMid = (2.0*cubicArcLength - lvRadius*radiansPerElementAroundLVFreeWall)/(elementsCountAroundVSeptum - 1)
        #elementLengthEnd = 0.5*(lvRadius*radiansPerElementAroundLVFreeWall + elementLengthMid)
        length = (elementsCountAroundVSeptum % 2)*elementLengthMid*0.5
        x = []
        dx_ds1 = []
        for n1 in range(elementsCountAroundVSeptum//2):
            xi = length/cubicArcLength
            pos1 = interp.interpolateCubicHermite(v1, d1, v2, d2, xi)
            x.append([ v for v in pos1 ])
            deriv1 = interp.interpolateCubicHermiteDerivative(v1, d1, v2, d2, xi)
            scale = elementLengthMid/vector.magnitude(deriv1)
            dx_ds1.append([ d*scale for d in deriv1 ])
            length += elementLengthMid
        return x, dx_ds1


def getRvOuterPoints(rvArcAroundRadians, lvRadius, rvAddWidthRadius, rvAddCrossRadius, elementsCountAroundRV, dEndMag, z, xi, mirror=False):
    '''
    Get array of points and derivatives around half of RV from +x axis anticlockwise.
    LV is assumed to be circular, centred at x, y = (0,0)
    :param rvArcAroundRadians: Angle in radians around outside of RV to tangent nodes where it meets LV.
    :param lvRadius: Radius of the LV.
    :param rvAddWidthRadius: Radius to add to LV to get RV width at centre.
    :param rvAddCrossRadius: Additional radius to add only laterally around septum.
    :param elementsCountAroundRV: Number of elements around RV.
    :param dEndMag: Magnitude of the external derivative at LV junction.
    :param z: Z coordinate to give to all values of x[]. dx_ds1[2] is all zero.
    :param xi: xi value ranging from 0 at bottom of RV to 1 at base.
    :param mirror: Set to True to mirror points, otherwise only second half of RV returned.
    :return: Arrays x[], dx_ds1[].
    '''
    elementsCountEnd = 0 if (xi <= 0.0) else 1
    #print('\ngetRvOuterPoints', rvArcAroundRadians, lvRadius, rvAddWidthRadius, rvAddCrossRadius, elementsCountAroundRV, dEndMag, z)
    startRadians = 0.5*rvArcAroundRadians
    ellipseEndRadians = startRadians - xi*0.5*rvArcAroundRadians/elementsCountAroundRV
    b = lvRadius + rvAddCrossRadius  # cross radius
    if rvAddCrossRadius >= 0.0:
        theta = math.pi - math.asin(lvRadius*math.sin(ellipseEndRadians)/b)
    else:
        theta = math.pi/2.0
    #print('phi',ellipseEndRadians)
    #print('theta',theta)
    a = (lvRadius*(math.cos(ellipseEndRadians) - 1.0) - rvAddWidthRadius)/(math.cos(theta) - 1.0)  # width radius
    cx = lvRadius + rvAddWidthRadius - a
    #print(ellipseEndRadians,'a',a,'b',b,'cx',cx)
    # get cubic curve joining half ellipse with to end point, first with unit derivatives then compute arc length
    x1 = ( cx, b )
    d1 = ( -1.0, 0.0 )
    x2 = ( lvRadius*math.cos(startRadians), lvRadius*math.sin(startRadians) )
    d2 = ( -math.sin(startRadians), math.cos(startRadians) )
    cubicArcLength = interp.computeCubicHermiteArcLength(x1, d1, x2, d2, True)
    d1 = ( d1[0]*cubicArcLength, d1[1]*cubicArcLength )
    d2 = ( d2[0]*cubicArcLength, d2[1]*cubicArcLength )
    quarterEllipsePerimeter = 0.25*getApproximateEllipsePerimeter(a, b)
    halfLength = cubicArcLength + quarterEllipsePerimeter
    #print('getRvOuterPoints. length cubic', cubicArcLength, ', quarterEllipsePerimeter', quarterEllipsePerimeter, ', halfLength',halfLength)

    if elementsCountEnd == 0:
        elementLengthTrans = elementLengthMid = (2.0*halfLength - dEndMag)/(elementsCountAroundRV - 1)
        elementLengthEnd = 0.5*(elementLengthMid + dEndMag)
    else:
        elementLengthEnd = dEndMag
        elementLengthMid = (2.0*halfLength - elementLengthEnd*3.0)/(elementsCountAroundRV - 3)
        elementLengthTrans = 0.5*(elementLengthEnd + elementLengthMid)
    #print('elementLengthMid',elementLengthMid,', elementLengthEnd',elementLengthEnd)
    length = 0.0
    angle = 0.0
    if (elementsCountAroundRV % 2) == 1:
        length += 0.5*elementLengthMid
        angle = updateEllipseAngleByArcLength(a, b, angle, length)
    x = []
    dx_ds1 = []
    nLimit = elementsCountAroundRV//2
    for n in range(nLimit):
        dMag = elementLengthMid
        if n < (nLimit - 2):
            elementLength = elementLengthMid
        elif n == (nLimit - 2):
            elementLength = elementLengthTrans
        else:
            elementLength = elementLengthEnd
            if elementsCountEnd > 0:
                dMag = elementLengthEnd
        if length <= quarterEllipsePerimeter:
            cosAngle = math.cos(angle)
            sinAngle = math.sin(angle)
            v = ( cx + a*cosAngle, b*sinAngle )
            t = ( -a*sinAngle, b*cosAngle )
            angle = updateEllipseAngleByArcLength(a, b, angle, elementLength)
        else:
            xi = (length - quarterEllipsePerimeter)/cubicArcLength
            v = interp.interpolateCubicHermite(x1, d1, x2, d2, xi)
            t = interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, xi)
        x.append([ v[0], v[1], z ])
        scale = dMag/math.sqrt(t[0]*t[0] + t[1]*t[1])
        dx_ds1.append([ t[0]*scale, t[1]*scale, 0.0 ])
        length += elementLength
    #print('getRvOuterPoints. length',length,', overshoot', length - halfLength)
    if mirror:
        xm = []
        dxm_ds1 = []
        for n in range(elementsCountAroundRV//2 - 1, -(elementsCountAroundRV%2), -1):
            xm.append([ x[n][0], -x[n][1], x[n][2] ])
            dxm_ds1.append([ -dx_ds1[n][0], dx_ds1[n][1], dx_ds1[n][2] ])
        x = xm + x
        dx_ds1 = dxm_ds1 + dx_ds1
    return x, dx_ds1


def getRVOuterSize(xiUpWidth, xiUpCross, rvWidthRadius, rvExtraCrossRadiusBase):
    if xiUpWidth < 0.0:
        #print('getRVOuterSize NEGATIVE xiUpWidth', xiUpWidth)
        return 0.0, 0.0
    if xiUpCross < 0.0:
        xiUpCross = 0.0
    xiUpFast = 1.0 - (1.0 - xiUpWidth)*(1.0 - xiUpWidth)
    xiUpSlow = xiUpCross
    rvAddWidthRadius = interp.interpolateCubicHermite([0.0], [0.0], [rvWidthRadius], [0.0], xiUpFast)[0]
    rvAddCrossRadius = interp.interpolateCubicHermite([0.0], [0.0], [rvExtraCrossRadiusBase], [0.0], xiUpSlow)[0]
    #print('getRVOuterSize xiUpWidth', xiUpWidth, ', addWidth', rvAddWidthRadius, ', addCross', rvAddCrossRadius)
    return rvAddWidthRadius, rvAddCrossRadius
