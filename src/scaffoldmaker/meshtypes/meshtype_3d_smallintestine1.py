"""
Generates a 3-D small intestine mesh along the central line,
with variable numbers of elements around, along and through
wall, with variable radius and thickness along.
"""

import copy
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.annotation.smallintestine_terms import get_smallintestine_term
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils.tubemesh import CylindricalSegmentTubeMeshInnerPoints
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from opencmiss.zinc.node import Node

class MeshType_3d_smallintestine1(Scaffold_base):
    '''
    Generates a 3-D small intestine mesh with variable numbers
    of elements around, along the central line, and through wall.
    The small intestine is created by a function that generates
    a small intestine segment and uses tubemesh to map the segment
    along a central line profile.
    '''

    centralPathDefaultScaffoldPackages = {
        'Mouse 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 45
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2  ], [
                [ [  -2.3, 18.5,  -4.4 ], [ -4.2, -0.8,   3.7 ], [  0.0,  5.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -8.6, 16.3,  -0.4 ], [ -7.1, -2.7,   1.6 ], [  0.0,  5.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -18.3, 12.6,  -1.5 ], [ -6.4, -1.7,  -3.8 ], [  0.0,  5.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -15.6, 13.7,  -6.1 ], [  7.0,  2.1,  -1.8 ], [  0.0,  5.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -9.3, 14.8,  -4.9 ], [  4.7,  0.7,   1.8 ], [  0.0,  5.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -3.9, 15.7,  -3.0 ], [  4.3,  0.7,   2.0 ], [  0.0,  5.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -3.4, 13.4,  -2.8 ], [ -4.1, -0.7,  -1.7 ], [  0.6, -2.0,  0.3 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -7.6, 12.4,  -4.6 ], [ -3.7, -0.8,  -0.9 ], [  0.0, -2.1,  0.1 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -11.6, 11.6,  -5.7 ], [ -4.2, -0.7,  -0.2 ], [  0.0, -1.9,  0.1 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -16.5, 11.7,  -3.9 ], [ -1.0,  0.2,   5.8 ], [  0.3, -1.4, -0.2 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -12.5, 11.7,  -1.4 ], [  3.6,  0.1,   0.6 ], [ -0.1, -1.4, -0.7 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -6.8, 11.8,  -0.6 ], [  2.9,  0.0,   0.7 ], [ -0.7, -1.2, -0.9 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -6.4,  9.8,  -1.6 ], [ -2.9, -0.3,  -1.4 ], [ -0.9,  1.6,  0.4 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -9.5,  9.5,  -2.9 ], [ -4.6,  0.0,  -1.8 ], [ -0.5,  1.7,  0.7 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -14.3,  9.4,  -4.6 ], [ -3.4,  0.1,  -1.6 ], [ -0.1,  1.6,  0.5 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -19.0,  9.4,  -2.9 ], [  0.3,  0.2,   6.7 ], [  0.0,  1.8,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -14.5,  9.7,   0.2 ], [  3.6, -1.2,   1.0 ], [  0.7,  2.2,  0.9 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -12.6,  7.7,   0.7 ], [  0.6, -2.7,   0.2 ], [  1.7,  0.8,  1.6 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -13.1,  3.8,   0.3 ], [ -4.0, -3.6,  -1.5 ], [  1.0, -1.1,  1.5 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -15.2,  5.1,  -0.8 ], [  6.0,  6.9,   1.8 ], [ -0.9, -0.4,  0.3 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -17.3,  6.9,  -1.0 ], [ -2.5,  0.0,  -0.4 ], [ -1.6,  0.0,  1.7 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -19.9,  6.8,  -2.5 ], [ -1.5, -1.1,  -3.4 ], [ -1.8, -0.7,  0.4 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -17.2,  6.3,  -5.1 ], [  4.0,  0.8,  -1.3 ], [ -0.4, -1.3, -1.5 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -12.2,  7.8,  -6.8 ], [  4.8,  1.7,  -0.3 ], [  0.1, -0.6, -1.6 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -7.9,  9.6,  -6.5 ], [  3.7,  1.7,   0.7 ], [  0.5, -0.6, -1.6 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -3.8, 10.3,  -5.5 ], [  3.8, -2.7,  -0.1 ], [ -1.3, -0.3, -3.1 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -5.3,  7.6,  -6.4 ], [ -3.5, -1.0,  -1.3 ], [ -0.4,  1.1, -1.5 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -9.0,  6.4,  -7.3 ], [ -3.2, -1.3,   1.9 ], [ -0.9,  0.5, -1.1 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -11.6,  4.0,  -2.0 ], [  5.6, -0.2,   4.3 ], [ -1.8, -0.3,  1.2 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -5.9,  5.0,  -3.1 ], [  4.1,  1.2,  -1.6 ], [  1.0, -0.3,  1.9 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -2.5,  6.0,  -3.8 ], [  3.6,  0.7,   3.2 ], [ -1.7, -1.2,  1.2 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -4.1,  3.2,  -0.4 ], [ -3.5, -1.7,   2.6 ], [ -1.0,  0.4, -1.1 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -9.7,  1.7,   2.3 ], [ -7.9, -1.0,   1.0 ], [ -0.6,  1.0, -1.6 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -19.0,  0.6,  -0.4 ], [  0.2,  3.7,  -6.8 ], [  2.6,  0.8,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -13.9,  2.3,  -5.8 ], [  4.4,  0.6,  -1.1 ], [  0.5, -0.1,  2.4 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -7.7,  1.2,  -4.6 ], [  3.9, -3.4,   1.5 ], [ -0.5,  0.0,  2.2 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -4.8, -4.0,  -1.3 ], [ -4.2, -3.3,   3.1 ], [ -1.5,  2.5,  2.3 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -10.9, -6.1,  -0.6 ], [ -5.3, -1.2,  -0.9 ], [ -0.5,  1.9,  1.3 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -19.9, -6.4,  -5.5 ], [ -0.3,  1.7, -10.5 ], [ -0.9,  2.0,  0.7 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -10.7, -3.2,  -8.8 ], [  7.8,  0.4,   0.1 ], [ -0.2,  2.1,  0.2 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -1.2, -1.9,  -7.3 ], [  0.8,  8.1,   2.5 ], [ -3.8,  0.1, -0.6 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -6.3,  0.5,  -8.1 ], [ -9.8, -1.2,   0.5 ], [  0.2, -1.9, -0.5 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -16.0, -0.7,  -7.4 ], [ -7.6,  1.2,   1.5 ], [ -0.1, -2.0, -0.9 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -20.5,  2.3,  -6.1 ], [  3.5,  7.2,  -2.9 ], [ -2.3,  0.2, -0.2 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -11.4,  2.6, -10.1 ], [ 10.4,  1.5,  -0.2 ], [ -0.2,  1.3, -2.2 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -3.8,  4.2,  -7.3 ], [  3.5,  0.9,   2.7 ], [  0.2,  1.0, -2.6 ], [ 0.0, 0.0, 0.5 ] ] ] )
            } )
        }

    @staticmethod
    def getName():
        return '3D Small Intestine 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Mouse 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 1']
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of segments': 100,
            'Number of elements around': 8,
            'Number of elements along segment': 4,
            'Number of elements through wall': 1,
            'Duodenum length': 25.0,
            'Jejunum length': 240.0,
            'Ileum length': 10.0,
            'Duodenum inner radius': 0.6,
            'Duodenum-jejunum inner radius': 1.0,
            'Jejunum-ileum inner radius': 1.0,
            'Ileum inner radius': 1.0,
            'Wall thickness': 0.1,
            'Use cross derivatives': False,
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements along': 1,
            'Refine number of elements through wall': 1
            }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of segments',
            'Number of elements around',
            'Number of elements along segment',
            'Number of elements through wall',
            'Duodenum length',
            'Jejunum length',
            'Ileum length',
            'Duodenum inner radius',
            'Duodenum-jejunum inner radius',
            'Jejunum-ileum inner radius',
            'Ileum inner radius',
            'Wall thickness',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall' ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [ MeshType_1d_path1 ]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
            'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        '''
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        '''
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + \
                ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Central path':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        for key in [
            'Number of segments',
            'Number of elements around',
            'Number of elements along segment',
            'Number of elements through wall',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Duodenum length',
            'Jejunum length',
            'Ileum length',
            'Duodenum inner radius',
            'Duodenum-jejunum inner radius',
            'Jejunum-ileum inner radius',
            'Ileum inner radius',
            'Wall thickness']:
            if options[key] < 0.0:
                options[key] = 0.0

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        """
        centralPath = options['Central path']
        segmentCount = options['Number of segments']
        elementsCountAround = options['Number of elements around']
        elementsCountAlongSegment = options['Number of elements along segment']
        elementsCountThroughWall = options['Number of elements through wall']
        duodenumLength = options['Duodenum length']
        jejunumLength = options['Jejunum length']
        duodenumInnerRadius = options['Duodenum inner radius']
        duodenumJejunumInnerRadius = options['Duodenum-jejunum inner radius']
        jejunumIleumInnerRadius = options['Jejunum-ileum inner radius']
        ileumInnerRadius = options['Ileum inner radius']
        wallThickness = options['Wall thickness']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])
        elementsCountAlong = int(elementsCountAlongSegment*segmentCount)
        startPhase = 0.0

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        # Central path
        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx, cd1, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)
        # for i in range(len(cx)):
        #     print(i, '[', cx[i], ',', cd1[i], ',', cd2[i],',', cd12[i], '],')
        del tmpRegion

        # find arclength of colon
        length = 0.0
        elementsCountIn = len(cx) - 1
        sd1 = interp.smoothCubicHermiteDerivativesLine(cx, cd1, fixAllDirections = True,
            magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
        for e in range(elementsCountIn):
            arcLength = interp.getCubicHermiteArcLength(cx[e], sd1[e], cx[e + 1], sd1[e + 1])
            # print(e+1, arcLength)
            length += arcLength
        segmentLength = length / segmentCount
        elementAlongLength = length / elementsCountAlong
        # print('Length = ', length)

        # Sample central path
        sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlongSegment*segmentCount)
        sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)

        # Generate variation of radius & tc width along length
        lengthList = [0.0, duodenumLength, duodenumLength + jejunumLength, length]
        innerRadiusList = [duodenumInnerRadius, duodenumJejunumInnerRadius, jejunumIleumInnerRadius, ileumInnerRadius]
        innerRadiusSegmentList, dInnerRadiusSegmentList = interp.sampleParameterAlongLine(lengthList, innerRadiusList,
                                                                                          segmentCount)

        # Create annotation groups for small intestine sections
        elementsAlongDuodenum = round(duodenumLength / elementAlongLength)
        elementsAlongJejunum = round(jejunumLength / elementAlongLength)
        elementsAlongIleum = elementsCountAlong - elementsAlongDuodenum - elementsAlongJejunum
        elementsCountAlongGroups = [elementsAlongDuodenum, elementsAlongJejunum, elementsAlongIleum]

        smallintestineGroup = AnnotationGroup(region, get_smallintestine_term("small intestine"))
        duodenumGroup = AnnotationGroup(region, get_smallintestine_term("duodenum"))
        jejunumGroup = AnnotationGroup(region, get_smallintestine_term("jejunum"))
        ileumGroup = AnnotationGroup(region, get_smallintestine_term("ileum"))

        annotationGroupAlong = [[smallintestineGroup, duodenumGroup],
                                [smallintestineGroup, jejunumGroup],
                                [smallintestineGroup, ileumGroup]]

        annotationGroupsAlong = []
        for i in range(len(elementsCountAlongGroups)):
            elementsCount = elementsCountAlongGroups[i]
            for n in range(elementsCount):
                annotationGroupsAlong.append(annotationGroupAlong[i])

        annotationGroupsAround = []
        for i in range(elementsCountAround):
            annotationGroupsAround.append([ ])

        annotationGroupsThroughWall = []
        for i in range(elementsCountThroughWall):
            annotationGroupsThroughWall.append([ ])

        xExtrude = []
        d1Extrude = []
        d2Extrude = []
        d3UnitExtrude = []

        # Create object
        smallIntestineSegmentTubeMeshInnerPoints = CylindricalSegmentTubeMeshInnerPoints(
            elementsCountAround, elementsCountAlongSegment, segmentLength,
            wallThickness, innerRadiusSegmentList, dInnerRadiusSegmentList, startPhase)

        for nSegment in range(segmentCount):
            # Create inner points
            xInner, d1Inner, d2Inner, transitElementList, segmentAxis, radiusAlongSegmentList = \
               smallIntestineSegmentTubeMeshInnerPoints.getCylindricalSegmentTubeMeshInnerPoints(nSegment)

            # Project reference point for warping onto central path
            start = nSegment*elementsCountAlongSegment
            end = (nSegment + 1)*elementsCountAlongSegment + 1
            sxRefList, sd1RefList, sd2ProjectedListRef, zRefList = \
                tubemesh.getPlaneProjectionOnCentralPath(xInner, elementsCountAround, elementsCountAlongSegment,
                                                         segmentLength, sx[start:end], sd1[start:end], sd2[start:end],
                                                         sd12[start:end])

            # Warp segment points
            xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList = tubemesh.warpSegmentPoints(
                xInner, d1Inner, d2Inner, segmentAxis, sxRefList, sd1RefList, sd2ProjectedListRef,
                elementsCountAround, elementsCountAlongSegment, zRefList, radiusAlongSegmentList,
                closedProximalEnd=False)

            # Store points along length
            xExtrude = xExtrude + (xWarpedList if nSegment == 0 else xWarpedList[elementsCountAround:])
            d1Extrude = d1Extrude + (d1WarpedList if nSegment == 0 else d1WarpedList[elementsCountAround:])

            # Smooth d2 for nodes between segments and recalculate d3
            if nSegment == 0:
                d2Extrude = d2Extrude + (d2WarpedList[:-elementsCountAround])
                d3UnitExtrude = d3UnitExtrude + (d3WarpedUnitList[:-elementsCountAround])
            else:
                xSecondFace = xWarpedList[elementsCountAround:elementsCountAround*2]
                d2SecondFace = d2WarpedList[elementsCountAround:elementsCountAround*2]
                for n1 in range(elementsCountAround):
                    nx = [xLastTwoFaces[n1], xLastTwoFaces[n1 + elementsCountAround], xSecondFace[n1]]
                    nd2 = [d2LastTwoFaces[n1], d2LastTwoFaces[n1 + elementsCountAround], d2SecondFace[n1]]
                    d2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True,
                                                                  fixEndDerivative = True)[1]
                    d2Extrude.append(d2)
                    d3Unit = vector.normalise(vector.crossproduct3(vector.normalise(d1LastTwoFaces[n1 + elementsCountAround]),
                                                                   vector.normalise(d2)))
                    d3UnitExtrude.append(d3Unit)
                d2Extrude = d2Extrude + \
                            (d2WarpedList[elementsCountAround:-elementsCountAround] if nSegment < segmentCount - 1 else
                             d2WarpedList[elementsCountAround:])
                d3UnitExtrude = d3UnitExtrude + \
                                (d3WarpedUnitList[elementsCountAround:-elementsCountAround] if nSegment < segmentCount - 1 else
                                 d3WarpedUnitList[elementsCountAround:])
            xLastTwoFaces = xWarpedList[-elementsCountAround*2:]
            d1LastTwoFaces = d1WarpedList[-elementsCountAround*2:]
            d2LastTwoFaces = d2WarpedList[-elementsCountAround*2:]

        # Create coordinates and derivatives
        xList, d1List, d2List, d3List, curvatureList = tubemesh.getCoordinatesFromInner(xExtrude, d1Extrude,
            d2Extrude, d3UnitExtrude, [wallThickness]*(elementsCountAlong+1),
            elementsCountAround, elementsCountAlong, elementsCountThroughWall, transitElementList)

        flatWidthList, xiList = smallIntestineSegmentTubeMeshInnerPoints.getFlatWidthAndXiList()

        # Create flat and texture coordinates
        xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture = tubemesh.createFlatAndTextureCoordinates(
            xiList, flatWidthList, length, wallThickness, elementsCountAround,
            elementsCountAlong, elementsCountThroughWall, transitElementList)

        # Create nodes and elements
        nextNodeIdentifier, nextElementIdentifier, annotationGroups = tubemesh.createNodesAndElements(
            region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture,
            elementsCountAround, elementsCountAlong, elementsCountThroughWall,
            annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
            firstNodeIdentifier, firstElementIdentifier, useCubicHermiteThroughWall, useCrossDerivatives,
            closedProximalEnd=False)

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong,
                                                       refineElementsCountThroughWall)
        return
