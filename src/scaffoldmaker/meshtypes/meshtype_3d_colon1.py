"""
Generates a 3-D colon mesh along the central line, with variable
numbers of elements around, along and through wall, with
variable radius and thickness along.
"""

import copy
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.annotation.colon_terms import get_colon_term
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import MeshType_3d_colonsegment1, ColonSegmentTubeMeshInnerPoints,\
    getTeniaColi, createFlatAndTextureCoordinatesTeniaColi, createNodesAndElementsTeniaColi
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from opencmiss.zinc.node import Node

class MeshType_3d_colon1(Scaffold_base):
    '''
    Generates a 3-D colon mesh with variable numbers
    of elements around, along the central line, and through wall.
    The colon is created by a function that generates a colon
    segment and uses tubemesh to map the segment along a central
    line profile.
    '''

    centralPathDefaultScaffoldPackages = {
        'Human 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 8
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2  ], [
                [ [   0.0,   0.0, 0.0 ], [ -50.7,  178.2, 0.0 ], [ -24.0,  -6.0, -12.0 ], [ -14.0,  -1.0, -12.0 ] ],
                [ [ -47.4, 188.6, 0.0 ], [ -19.3,  177.1, 0.0 ], [ -22.0,  -4.0,  -8.0 ], [  -4.0,  19.0,  22.0 ] ],
                [ [  -4.4, 396.5, 0.0 ], [ 206.0,   40.1, 0.0 ], [ -10.0,  20.0,   8.0 ], [  -6.0,   0.0,  51.0 ] ],
                [ [ 130.0, 384.1, 0.0 ], [ 130.8,  -40.5, 0.0 ], [  -5.0,   4.0,  29.0 ], [   0.0,   1.0,  24.0 ] ],
                [ [ 279.4, 383.0, 0.0 ], [ 118.0,   48.7, 0.0 ], [  -2.0,  10.0,  22.0 ], [   5.0,  25.0, -20.0 ] ],
                [ [ 443.9, 390.8, 0.0 ], [ 111.3,  -97.0, 0.0 ], [  10.0,  17.0,   6.0 ], [   1.0,  -6.0, -35.0 ] ],
                [ [ 475.2, 168.0, 0.0 ], [  -0.8, -112.4, 0.0 ], [  20.0,   0.0, -20.0 ], [  15.0,  -1.0, -10.0 ] ],
                [ [ 432.6, -32.3, 0.0 ], [ -90.5,  -59.0, 0.0 ], [   6.0,  -9.0, -14.0 ], [   8.0, -11.0, -13.0 ] ],
                [ [ 272.4,   7.5, 0.0 ], [ -79.0,   47.4, 0.0 ], [   1.0, -11.0, -18.0 ], [   4.0, -12.0, -12.0 ] ] ] )
            } ),
        'Human 2' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 8
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [   0.0,   0.0,    0.0 ], [ -34.7,  104.1, -18.1 ], [ -24.0,  -6.0, -12.0 ], [ -14.0,  -1.0, -12.0 ] ],
                [ [ -34.5,  114.0, -18.1 ], [   1.2,   86.6,  -3.4 ], [ -22.0,  -4.0,  -8.0 ], [  -4.0,  19.0,  22.0 ] ],
                [ [ -19.1,  218.5,   5.5 ], [  78.7,   -7.1,  94.5 ], [ -10.0,  20.0,   8.0 ], [  -6.0,   0.0,  51.0 ] ],
                [ [  82.5,  189.1,  94.2 ], [  84.5,    7.1,  71.6 ], [  -5.0,   4.0,  29.0 ], [   0.0,   1.0,  24.0 ] ],
                [ [ 226.6,  218.7,  85.7 ], [  95.0,   91.3, -58.5 ], [  -2.0,  10.0,  22.0 ], [   5.0,  25.0, -20.0 ] ],
                [ [ 325.5,  381.7, -57.9 ], [ 229.2,  -66.7, -20.4 ], [  10.0,  17.0,   6.0 ], [   1.0,  -6.0, -35.0 ] ],
                [ [ 354.0,  105.3, -24.4 ], [  -6.3, -143.7,  20.3 ], [  20.0,   0.0, -20.0 ], [  15.0,  -1.0, -10.0 ] ],
                [ [ 296.5, -121.2,  -0.6 ], [ -90.5,  -59.0,   0.0 ], [   6.0,  -9.0, -14.0 ], [   8.0, -11.0, -13.0 ] ],
                [ [ 169.8,  -73.4, -33.5 ], [ -72.2,   43.4, -27.4 ], [   1.0, -11.0, -18.0 ], [   4.0, -12.0, -12.0 ] ] ] )
            } ),
        'Mouse 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 7
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [   0.0,   0.0,  0.0 ], [  6.0, 12.0,  -2.0 ], [ 2.0,  1.0,  2.0 ], [ 6.0, 0.0, 3.0 ] ],
                [ [  -2.0,  11.0, -3.0 ], [ -8.0,  4.0,   9.0 ], [ 2.0,  2.0,  1.0 ], [ 0.0, 1.0, 2.0 ] ],
                [ [  -3.0,   2.0,  3.0 ], [ -4.0, -8.0,   0.0 ], [ 2.0, -1.0,  2.0 ], [ 1.0, 0.0, 2.0 ] ],
                [ [ -11.0,  -3.0, -4.0 ], [ -8.0, -3.0,  -7.0 ], [ 1.0, -2.0,  1.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -16.0,  -4.0,  0.0 ], [  4.0, -3.0,  14.0 ], [ 1.0, -3.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -7.0,  -8.0,  0.0 ], [  5.0, -1.0, -14.0 ], [ 0.0, -3.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -1.0,  -6.0, -1.0 ], [  2.0, -2.0,   9.0 ], [ 1.0, -3.0, -1.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -2.0, -14.0,  5.0 ], [ -2.0, -4.0,   2.0 ], [ 1.0, -2.0, -2.0 ], [ 0.0, 0.0, 0.5 ] ] ] )
            } ),
        'Mouse 2' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 4
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [   0.0,  0.0,   0.0 ], [  0.0,  0.0,  13.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   0.0,  0.0,  13.0 ], [  0.0,  2.0,  28.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -14.0, -2.0,  13.0 ], [  0.0, -3.0, -19.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -14.0, -1.0, -10.0 ], [  1.0,  1.0, -17.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -14.0,  0.0, -28.0 ], [  0.0,  0.0, -11.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ] ] )
            } ),
        'Pig 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 36
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [   -7.2,   83.3,  -20.7 ], [  -65.2,   -8.1,   7.6 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -68.5,   52.8,   -9.6 ], [  -40.1,  -36.1,  10.7 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -97.4,  -26.3,    5.7 ], [   18.0,  -93.2,  13.7 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -56.8,  -90.5,   14.1 ], [   65.5,  -41.4,   7.3 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   48.9, -100.8,   24.0 ], [  112.2,   40.1,  19.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  114.8,  -12.6,   38.7 ], [    8.2,   96.1,  14.2 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   60.3,   83.5,   43.7 ], [ -108.7,   54.1,  22.4 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -41.2,   90.7,   56.3 ], [  -89.0,  -32.4,  14.4 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -107.9,   -9.7,   76.6 ], [   11.1,  -94.4,  11.3 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -57.3,  -91.9,   81.3 ], [   71.2,  -31.2,   5.7 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   51.2,  -89.4,   97.2 ], [   99.1,   55.4,  12.9 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   91.6,    9.3,  103.6 ], [    4.7,   51.2,   3.4 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   61.6,  111.8,  109.6 ], [  -85.2,   46.1,   2.6 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -54.6,   91.9,  129.4 ], [  -92.7,  -55.0,  14.5 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -109.0,    5.6,  156.9 ], [   23.6, -108.2,  27.7 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -59.1,  -62.5,  170.8 ], [   74.0,  -20.1,  14.4 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   23.5,  -53.2,  179.7 ], [   84.6,   47.0,   6.9 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   62.3,   30.1,  187.5 ], [  -12.8,   58.0,   0.8 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   22.4,   45.2,  181.1 ], [  -23.6,  -34.5,  -7.4 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   -1.9,    4.9,  180.5 ], [  -41.3,  -30.9,   7.5 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -45.1,  -12.6,  194.4 ], [  -40.5,   -4.6,   6.9 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -71.7,   -2.2,  197.2 ], [  -25.2,   35.8,  -6.8 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -65.8,   42.1,  182.3 ], [   26.6,   37.6, -15.6 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -14.1,   81.2,  163.5 ], [   41.0,   10.3,  -9.5 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   61.7,   86.1,  156.4 ], [   77.9,  -40.7,   8.9 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   92.9,   20.5,  150.3 ], [    0.0,  -73.3,  -5.2 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   48.9,  -65.0,  142.8 ], [  -82.8,  -80.0,  -1.9 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -54.3,  -90.8,  134.0 ], [  -60.1,   26.4,  -8.2 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -89.9,   11.2,  115.0 ], [   34.9,  125.1, -27.9 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -17.4,   74.2,   91.1 ], [   78.8,   19.1, -15.4 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   43.4,   50.2,   73.3 ], [   30.2,  -36.0,  -9.9 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   62.4,   -5.1,   63.5 ], [   10.9,  -54.2,  -2.7 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   32.7,  -51.7,   56.1 ], [  -38.6,  -29.8,  -8.1 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -38.1,  -28.6,   46.8 ], [  -66.0,   83.3, -12.1 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   -5.4,   32.0,   26.0 ], [   48.1,   17.6, -21.4 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  146.6,  101.2,  -41.2 ], [   63.3,   35.3, -31.2 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  312.3,  199.5, -123.4 ], [   39.7,   24.3, -20.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ] ] )
            } ),
        'Pig 2': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 90.0,
                'Number of elements': 3
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                [ [  0.0, 0.0, 0.0 ], [ 30.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ],
                [ [ 30.0, 0.0, 0.0 ], [ 30.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ],
                [ [ 60.0, 0.0, 0.0 ], [ 30.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ],
                [ [ 90.0, 0.0, 0.0 ], [ 30.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ] ] )
            } ),
        }

    @staticmethod
    def getName():
        return '3D Colon 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Human 2',
            'Mouse 1',
            'Mouse 2',
            'Pig 1',
            'Pig 2']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Human 2' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 2']
        elif 'Mouse 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 1']
        elif 'Mouse 2' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 2']
        elif 'Pig 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Pig 1']
        elif 'Pig 2' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Pig 2']
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']
        if 'Mouse' in parameterSetName:
            segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName = 'Mouse 1')
        elif 'Pig' in parameterSetName:
            segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName = 'Pig 1')
        else:
            segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName = 'Human 1')
        options = {
            'Central path' : copy.deepcopy(centralPathOption),
            'Segment profile' : segmentProfileOption,
            'Number of segments': 30,
            'Start phase': 0.0,
            'Proximal length': 420.0,
            'Transverse length': 460.0,
            'Distal length': 620.0,
            'Proximal inner radius': 43.5,
            'Proximal tenia coli width': 10.0,
            'Proximal-transverse inner radius': 33.0,
            'Proximal-transverse tenia coli width': 10.0,
            'Transverse-distal inner radius': 29.0,
            'Transverse-distal tenia coli width': 10.0,
            'Distal inner radius': 31.5,
            'Distal tenia coli width': 10.0,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
            }
        if 'Human 2' in parameterSetName:
            options['Proximal length'] = 180.0
            options['Transverse length'] = 620.0
            options['Distal length'] = 700.0
        elif 'Mouse' in parameterSetName:
            options['Number of segments'] = 10
            options['Proximal length'] = 30.0
            options['Transverse length'] = 20.0
            options['Distal length'] = 25.0
            options['Proximal inner radius'] = 1.0
            options['Proximal tenia coli width'] = 0.5
            options['Proximal-transverse inner radius'] = 0.9
            options['Proximal-transverse tenia coli width'] = 0.7
            options['Transverse-distal inner radius'] = 0.7
            options['Transverse-distal tenia coli width'] = 1.0
            options['Distal inner radius'] = 0.7
            options['Distal tenia coli width'] = 1.0
        elif 'Pig 1' in parameterSetName:
            options['Number of segments'] = 120
            options['Proximal length'] = 3000.0
            options['Transverse length'] = 200.0
            options['Distal length'] = 200.0
            options['Proximal inner radius'] = 38.0
            options['Proximal tenia coli width'] = 5.0
            options['Proximal-transverse inner radius'] = 14.0
            options['Proximal-transverse tenia coli width'] = 4.0
            options['Transverse-distal inner radius'] = 10.5
            options['Transverse-distal tenia coli width'] = 3.0
            options['Distal inner radius'] = 8.0
            options['Distal tenia coli width'] = 1.5
        elif 'Pig 2' in parameterSetName:
            options['Number of segments'] = 3
            options['Proximal length'] = 30.0
            options['Transverse length'] = 30.0
            options['Distal length'] = 30.0
            options['Proximal inner radius'] = 16.0
            options['Proximal tenia coli width'] = 5.0
            options['Proximal-transverse inner radius'] = 16.0
            options['Proximal-transverse tenia coli width'] = 5.0
            options['Transverse-distal inner radius'] = 16.0
            options['Transverse-distal tenia coli width'] = 5.0
            options['Distal inner radius'] = 16.0
            options['Distal tenia coli width'] = 5.0
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Segment profile',
            'Number of segments',
            'Start phase',
            'Proximal length',
            'Transverse length',
            'Distal length',
            'Proximal inner radius',
            'Proximal tenia coli width',
            'Proximal-transverse inner radius',
            'Proximal-transverse tenia coli width',
            'Transverse-distal inner radius',
            'Transverse-distal tenia coli width',
            'Distal inner radius',
            'Distal tenia coli width',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall' ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [ MeshType_1d_path1 ]
        if optionName == 'Segment profile':
            return [ MeshType_3d_colonsegment1 ]
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
        if optionName == 'Segment profile':
            if not parameterSetName:
                parameterSetName = scaffoldType.getParameterSetNames()[0]
            return ScaffoldPackage(scaffoldType, defaultParameterSetName = parameterSetName)
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        if not options['Segment profile'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Segment profile'):
            options['Segment profile'] = cls.getOptionScaffoldPackage('Segment profile', MeshType_3d_colonsegment1)
        for key in [
            'Number of segments',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Proximal length',
            'Transverse length',
            'Distal length',
            'Proximal inner radius',
            'Proximal tenia coli width',
            'Proximal-transverse inner radius',
            'Proximal-transverse tenia coli width',
            'Transverse-distal inner radius',
            'Transverse-distal tenia coli width',
            'Distal inner radius',
            'Distal tenia coli width']:
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
        segmentProfile = options['Segment profile']
        segmentCount = options['Number of segments']
        startPhase = options['Start phase'] % 360.0
        proximalLength = options['Proximal length']
        transverseLength = options['Transverse length']
        proximalInnerRadius = options['Proximal inner radius']
        proximalTCWidth = options['Proximal tenia coli width']
        proximalTransverseInnerRadius = options['Proximal-transverse inner radius']
        proximalTransverseTCWidth = options['Proximal-transverse tenia coli width']
        transverseDistalInnerRadius = options['Transverse-distal inner radius']
        transverseDistalTCWidth = options['Transverse-distal tenia coli width']
        distalInnerRadius = options['Distal inner radius']
        distalTCWidth = options['Distal tenia coli width']
        segmentSettings = segmentProfile.getScaffoldSettings()

        elementsCountAroundTC = segmentSettings['Number of elements around tenia coli']
        elementsCountAroundHaustrum = segmentSettings['Number of elements around haustrum']
        cornerInnerRadiusFactor = segmentSettings['Corner inner radius factor']
        haustrumInnerRadiusFactor = segmentSettings['Haustrum inner radius factor']
        segmentLengthEndDerivativeFactor = segmentSettings['Segment length end derivative factor']
        segmentLengthMidDerivativeFactor = segmentSettings['Segment length mid derivative factor']
        tcCount = segmentSettings['Number of tenia coli']
        tcThickness = segmentSettings['Tenia coli thickness']
        elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum)*tcCount

        elementsCountAlongSegment = segmentSettings['Number of elements along segment']
        elementsCountThroughWall = segmentSettings['Number of elements through wall']
        wallThickness = segmentSettings['Wall thickness']
        useCrossDerivatives = segmentSettings['Use cross derivatives']
        useCubicHermiteThroughWall = not(segmentSettings['Use linear through wall'])
        elementsCountAlong = int(elementsCountAlongSegment*segmentCount)

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        # Central path
        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx, cd1, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)
        # for i in range(len(cx)):
        #     print(i, '[', cx[i], ',', cd1[i], ',', cd2[i], ',', cd12[i], '],')
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
        # print('Length = ', length)
        elementAlongLength = length / elementsCountAlong

        # Sample central path
        sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlong)
        sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)

        # Generate variation of radius & tc width along length
        lengthList = [0.0, proximalLength, proximalLength + transverseLength, length]
        innerRadiusList = [proximalInnerRadius, proximalTransverseInnerRadius,
                           transverseDistalInnerRadius, distalInnerRadius]
        innerRadiusAlongElementList, dInnerRadiusAlongElementList = interp.sampleParameterAlongLine(lengthList,
                                                                                                    innerRadiusList,
                                                                                                    elementsCountAlong)

        tcWidthList = [proximalTCWidth, proximalTransverseTCWidth, transverseDistalTCWidth, distalTCWidth]
        tcWidthAlongElementList, dTCWidthAlongElementList = interp.sampleParameterAlongLine(lengthList,
                                                                                            tcWidthList,
                                                                                            elementsCountAlong)

        # Account for reduced haustrum appearance in transverse and distal pig colon
        if tcCount == 2:
            haustrumInnerRadiusFactorList = [haustrumInnerRadiusFactor, haustrumInnerRadiusFactor*0.75,
                                             haustrumInnerRadiusFactor*0.5, haustrumInnerRadiusFactor*0.2]
            haustrumInnerRadiusFactorAlongElementList = \
                interp.sampleParameterAlongLine(lengthList, haustrumInnerRadiusFactorList, elementsCountAlong)[0]
        else:
            haustrumInnerRadiusFactorAlongElementList = [haustrumInnerRadiusFactor]*(elementsCountAlong+1)

        # Create annotation groups for colon sections
        elementsAlongInProximal = round(proximalLength/elementAlongLength)
        elementsAlongInTransverse = round(transverseLength/elementAlongLength)
        elementsAlongInDistal = elementsCountAlong - elementsAlongInProximal - elementsAlongInTransverse
        elementsCountAlongGroups = [elementsAlongInProximal, elementsAlongInTransverse, elementsAlongInDistal]

        colonGroup = AnnotationGroup(region, get_colon_term("colon"))

        if tcCount == 1:
            proximalGroup = AnnotationGroup(region, get_colon_term("proximal colon"))
            transverseGroup = AnnotationGroup(region, get_colon_term("transverse colon"))
            distalGroup = AnnotationGroup(region, get_colon_term("distal colon"))
            annotationGroupAlong = [[colonGroup, proximalGroup],
                                    [colonGroup, transverseGroup],
                                    [colonGroup, distalGroup]]

        elif tcCount == 2:
            spiralGroup = AnnotationGroup(region, get_colon_term("spiral colon"))
            transverseGroup = AnnotationGroup(region, get_colon_term("transverse colon"))
            distalGroup = AnnotationGroup(region, get_colon_term("distal colon"))
            annotationGroupAlong = [[colonGroup, spiralGroup],
                                    [colonGroup, transverseGroup],
                                    [colonGroup, distalGroup]]

        elif tcCount == 3:
            ascendingGroup = AnnotationGroup(region, get_colon_term("ascending colon"))
            transverseGroup = AnnotationGroup(region, get_colon_term("transverse colon"))
            descendingGroup = AnnotationGroup(region, get_colon_term("descending colon"))
            annotationGroupAlong = [[colonGroup, ascendingGroup],
                                    [colonGroup, transverseGroup],
                                    [colonGroup, descendingGroup]]

        annotationGroupsAlong = []
        for i in range(len(elementsCountAlongGroups)):
            elementsCount = elementsCountAlongGroups[i]
            for n in range(elementsCount):
                annotationGroupsAlong.append(annotationGroupAlong[i])

        annotationGroupsThroughWall = []
        for i in range(elementsCountThroughWall):
            annotationGroupsThroughWall.append([ ])

        xExtrude = []
        d1Extrude = []
        d2Extrude = []
        d3UnitExtrude = []
        sxRefExtrudeList = []

        # Create object
        colonSegmentTubeMeshInnerPoints = ColonSegmentTubeMeshInnerPoints(
            region, elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongSegment,
            tcCount, segmentLengthEndDerivativeFactor, segmentLengthMidDerivativeFactor,
            segmentLength, wallThickness, cornerInnerRadiusFactor, haustrumInnerRadiusFactorAlongElementList,
            innerRadiusAlongElementList, dInnerRadiusAlongElementList, tcWidthAlongElementList,
            startPhase)

        for nSegment in range(segmentCount):
            # Create inner points
            xInner, d1Inner, d2Inner, transitElementList, segmentAxis, annotationGroupsAround \
                = colonSegmentTubeMeshInnerPoints.getColonSegmentTubeMeshInnerPoints(nSegment)

            # Project reference point for warping onto central path
            start = nSegment * elementsCountAlongSegment
            end = (nSegment + 1) * elementsCountAlongSegment + 1
            sxRefList, sd1RefList, sd2ProjectedListRef, zRefList = \
                tubemesh.getPlaneProjectionOnCentralPath(xInner, elementsCountAround, elementsCountAlongSegment,
                                                         segmentLength, sx[start:end], sd1[start:end], sd2[start:end],
                                                         sd12[start:end])

            # Warp segment points
            xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList = tubemesh.warpSegmentPoints(
                xInner, d1Inner, d2Inner, segmentAxis, sxRefList, sd1RefList, sd2ProjectedListRef,
                elementsCountAround, elementsCountAlongSegment, zRefList, innerRadiusAlongElementList[start:end],
                closedProximalEnd=False)

            # Store points along length
            xExtrude +=  xWarpedList if nSegment == 0 else xWarpedList[elementsCountAround:]
            d1Extrude += d1WarpedList if nSegment == 0 else d1WarpedList[elementsCountAround:]
            d2Extrude += d2WarpedList if nSegment == 0 else d2WarpedList[elementsCountAround:]
            d3UnitExtrude += d3WarpedUnitList if nSegment == 0 else d3WarpedUnitList[elementsCountAround:]
            sxRefExtrudeList += sxRefList if nSegment == 0 else sxRefList[elementsCountAround:]

        contractedWallThicknessList = colonSegmentTubeMeshInnerPoints.getContractedWallThicknessList()

        # Create coordinates and derivatives
        xList, d1List, d2List, d3List, curvatureList = tubemesh.getCoordinatesFromInner(xExtrude, d1Extrude,
            d2Extrude, d3UnitExtrude, contractedWallThicknessList,
            elementsCountAround, elementsCountAlong, elementsCountThroughWall, transitElementList)

        relaxedLengthList, xiList = colonSegmentTubeMeshInnerPoints.getRelaxedLengthAndXiList()

        closedProximalEnd = False

        if tcThickness > 0:
            tubeTCWidthList = colonSegmentTubeMeshInnerPoints.getTubeTCWidthList()
            xList, d1List, d2List, d3List, annotationArrayAround = getTeniaColi(
                region, xList, d1List, d2List, d3List, curvatureList, tcCount, elementsCountAroundTC,
                elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
                tubeTCWidthList, tcThickness, sxRefExtrudeList, annotationGroupsAround,
                closedProximalEnd)

            # Create flat and texture coordinates
            xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture = createFlatAndTextureCoordinatesTeniaColi(
                xiList, relaxedLengthList, length, wallThickness, tcCount, tcThickness,
                elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong,
                elementsCountThroughWall, transitElementList, closedProximalEnd)

            # Create nodes and elements
            nextNodeIdentifier, nextElementIdentifier, annotationGroups = createNodesAndElementsTeniaColi(
                region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture,
                elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
                tcCount, annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
                firstNodeIdentifier, firstElementIdentifier, useCubicHermiteThroughWall, useCrossDerivatives,
                closedProximalEnd)

        else:
            # Create flat and texture coordinates
            xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture = tubemesh.createFlatAndTextureCoordinates(
                xiList, relaxedLengthList, length, wallThickness, elementsCountAround,
                elementsCountAlong, elementsCountThroughWall, transitElementList)

            # Create nodes and elements
            nextNodeIdentifier, nextElementIdentifier, annotationGroups = tubemesh.createNodesAndElements(
                region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture,
                elementsCountAround, elementsCountAlong, elementsCountThroughWall,
                annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
                firstNodeIdentifier, firstElementIdentifier, useCubicHermiteThroughWall, useCrossDerivatives,
                closedProximalEnd)

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
