"""
Generates a 3-D stomach mesh along the central line, with variable
numbers of elements around esophagus and duodenum, along and through
wall, with variable radius and thickness along.
"""

from __future__ import division
import math
import copy
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, mergeAnnotationGroups, \
    getAnnotationGroupForTerm, findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.stomach_terms import get_stomach_term
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldStoredString, findOrCreateFieldStoredMeshLocation, findOrCreateFieldNodeGroup
from opencmiss.utils.zinc.finiteelement import get_element_node_identifiers
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_ostium1 import MeshType_3d_ostium1, generateOstiumMesh
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.utils.bifurcation import get_bifurcation_triple_point
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eft_utils import scaleEftNodeValueLabels, setEftScaleFactorIds, remapEftNodeValueLabel, \
    remapEftNodeValueLabelsVersion
from scaffoldmaker.utils.geometry import createEllipsePoints
from scaffoldmaker.utils.tracksurface import TrackSurface
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import vector


class MeshType_3d_stomach1(Scaffold_base):
    """
    Generates a 3-D stomach mesh with variable numbers of elements around the esophagus and duodenum,
    along the central line, and through wall. The stomach is created using a central path as the longitudinal axis
    of the stomach. D2 of the central path points to the greater curvature of the stomach and magnitude of D2 and D3
    are the radii of the stomach in the respective direction. D2 on the first node of the central path provide
    the radius of the fundus dome.
    """
    centralPathDefaultScaffoldPackages = {
        'Human 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 7
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [ [  70.8, 72.3, 0.0 ], [ -17.0, -47.9, 0.0 ], [  39.3, -13.9, 0.0 ], [  -3.7, -5.9, 0.0 ], [ 0.0, 0.0, 38.7 ], [ 0.0, 0.0,  4.5] ],
                [ [  51.1, 30.7, 0.0 ], [ -22.3, -34.9, 0.0 ], [  32.2, -20.5, 0.0 ], [ -10.6, -7.2, 0.0 ], [ 0.0, 0.0, 41.4 ], [ 0.0, 0.0,  0.9] ],
                [ [  27.8,  2.7, 0.0 ], [ -27.7, -18.4, 0.0 ], [  18.7, -28.2, 0.0 ], [ -15.1, -2.0, 0.0 ], [ 0.0, 0.0, 40.9 ], [ 0.0, 0.0, -5.0] ],
                [ [  -0.3, -5.8, 0.0 ], [ -28.0,  -2.6, 0.0 ], [   2.3, -25.6, 0.0 ], [ -12.9,  4.6, 0.0 ], [ 0.0, 0.0, 32.3 ], [ 0.0, 0.0, -9.4] ],
                [ [ -26.5, -2.9, 0.0 ], [ -20.6,   8.0, 0.0 ], [  -7.5, -19.2, 0.0 ], [  -6.4,  8.0, 0.0 ], [ 0.0, 0.0, 22.1 ], [ 0.0, 0.0, -6.8] ],
                [ [ -40.6,  7.4, 0.0 ], [ -11.0,  12.5, 0.0 ], [ -11.6, -10.2, 0.0 ], [  -0.4,  7.8, 0.0 ], [ 0.0, 0.0, 17.6 ], [ 0.0, 0.0, -5.9] ],
                [ [ -48.1, 21.0, 0.0 ], [  -6.4,  15.7, 0.0 ], [  -8.6,  -3.5, 0.0 ], [   0.4,  4.4, 0.0 ], [ 0.0, 0.0, 10.5 ], [ 0.0, 0.0, -3.1] ],
                [ [ -52.9, 38.6, 0.0 ], [  -3.2,  19.3, 0.0 ], [ -11.1,  -1.9, 0.0 ], [  -5.5, -1.1, 0.0 ], [ 0.0, 0.0, 12.0 ], [ 0.0, 0.0,  6.1] ] ] ),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '2-3',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '4-5',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6',
                    'name': get_stomach_term('pylorus')[0],
                    'ontId': get_stomach_term('pylorus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7',
                    'name': get_stomach_term('duodenum')[0],
                    'ontId': get_stomach_term('duodenum')[1]
                }]
            }),
        'Mouse 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 8
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [ [  0.9,  3.7, 0.0 ], [ -0.8, -3.6, 0.0 ], [  3.2, -0.6, 0.0 ], [ -1.3, -0.5, 0.0 ], [ 0.0, 0.0, 2.6 ], [ 0.0, 0.0,  0.9 ] ], 
                [ [ -0.1,  0.7, 0.0 ], [ -1.2, -2.4, 0.0 ], [  2.0, -1.5, 0.0 ], [ -1.1, -1.3, 0.0 ], [ 0.0, 0.0, 3.1 ], [ 0.0, 0.0,  0.1 ] ], 
                [ [ -1.4, -1.1, 0.0 ], [ -1.6, -1.1, 0.0 ], [  1.0, -3.0, 0.0 ], [ -1.3, -0.8, 0.0 ], [ 0.0, 0.0, 3.0 ], [ 0.0, 0.0, -0.2 ] ], 
                [ [ -2.9, -1.6, 0.0 ], [ -1.6,  0.2, 0.0 ], [ -0.6, -3.3, 0.0 ], [ -1.4,  0.2, 0.0 ], [ 0.0, 0.0, 2.8 ], [ 0.0, 0.0, -0.1 ] ], 
                [ [ -4.3, -0.8, 0.0 ], [ -1.2,  1.1, 0.0 ], [ -1.8, -2.5, 0.0 ], [ -0.8,  1.1, 0.0 ], [ 0.0, 0.0, 2.9 ], [ 0.0, 0.0, -0.1 ] ], 
                [ [ -5.2,  0.6, 0.0 ], [ -0.8,  1.6, 0.0 ], [ -2.2, -1.1, 0.0 ], [  0.2,  1.1, 0.0 ], [ 0.0, 0.0, 2.5 ], [ 0.0, 0.0, -0.7 ] ], 
                [ [ -5.9,  2.3, 0.0 ], [ -0.5,  1.3, 0.0 ], [ -1.3, -0.4, 0.0 ], [  0.6,  0.3, 0.0 ], [ 0.0, 0.0, 1.4 ], [ 0.0, 0.0, -0.7 ] ], 
                [ [ -6.2,  3.2, 0.0 ], [ -0.4,  0.9, 0.0 ], [ -0.8, -0.3, 0.0 ], [  0.1, -0.0, 0.0 ], [ 0.0, 0.0, 0.9 ], [ 0.0, 0.0, -0.2 ] ], 
                [ [ -6.8,  4.1, 0.0 ], [ -0.7,  0.9, 0.0 ], [ -1.1, -0.5, 0.0 ], [ -0.7, -0.4, 0.0 ], [ 0.0, 0.0, 1.1 ], [ 0.0, 0.0,  0.6 ] ] ] ),

                'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-3',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '4-5',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7',
                    'name': get_stomach_term('pylorus')[0],
                    'ontId': get_stomach_term('pylorus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '8',
                    'name': get_stomach_term('duodenum')[0],
                    'ontId': get_stomach_term('duodenum')[1]
                }]
        }),
        'Rat 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 8
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [ [  12.4, 11.9, 0.0 ], [ -0.1, -14.6, 0.0 ], [  8.7, -1.4, 0.0 ], [ -0.7, -3.5, 0.0 ], [ 0.0, 0.0, 7.4 ], [ 0.0, 0.0,  1.5] ],
                [ [   9.5, -0.2, 0.0 ], [ -5.7,  -9.0, 0.0 ], [  6.6, -5.1, 0.0 ], [ -3.5, -3.9, 0.0 ], [ 0.0, 0.0, 8.5 ], [ 0.0, 0.0,  0.7] ],
                [ [   2.8, -5.5, 0.0 ], [ -6.6,  -2.6, 0.0 ], [  2.1, -9.1, 0.0 ], [ -4.0, -1.6, 0.0 ], [ 0.0, 0.0, 9.0 ], [ 0.0, 0.0,  0.1] ],
                [ [  -2.6, -6.0, 0.0 ], [ -5.7,   0.7, 0.0 ], [ -1.6, -9.2, 0.0 ], [ -4.0,  0.5, 0.0 ], [ 0.0, 0.0, 8.8 ], [ 0.0, 0.0, -0.4] ],
                [ [  -8.2, -4.1, 0.0 ], [ -5.3,   3.9, 0.0 ], [ -6.0, -8.0, 0.0 ], [ -2.6,  2.7, 0.0 ], [ 0.0, 0.0, 8.1 ], [ 0.0, 0.0, -1.2] ],
                [ [ -12.3,  1.5, 0.0 ], [ -1.5,   7.0, 0.0 ], [ -6.4, -3.5, 0.0 ], [  1.5,  4.0, 0.0 ], [ 0.0, 0.0, 6.2 ], [ 0.0, 0.0, -2.8] ],
                [ [ -10.8,  8.5, 0.0 ], [  0.6,   5.4, 0.0 ], [ -2.9,  0.0, 0.0 ], [  0.8,  1.0, 0.0 ], [ 0.0, 0.0, 2.4 ], [ 0.0, 0.0, -0.6] ],
                [ [ -10.7, 12.2, 0.0 ], [ -0.3,   3.2, 0.0 ], [ -3.5, -0.3, 0.0 ], [ -0.3, -0.1, 0.0 ], [ 0.0, 0.0, 3.4 ], [ 0.0, 0.0,  0.4] ],
                [ [ -11.3, 14.8, 0.0 ], [ -0.9,   2.0, 0.0 ], [ -3.5, -0.3, 0.0 ], [  0.3,  0.1, 0.0 ], [ 0.0, 0.0, 3.4 ], [ 0.0, 0.0, -0.4] ] ] ),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-3',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '4-5',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7',
                    'name': get_stomach_term('pylorus')[0],
                    'ontId': get_stomach_term('pylorus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '8',
                    'name': get_stomach_term('duodenum')[0],
                    'ontId': get_stomach_term('duodenum')[1]
                }]
        }),
    }

    ostiumDefaultScaffoldPackages = {
        'Human 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 1,  # not implemented for > 1
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 25.0,
                'Ostium length': 15.0,
                'Ostium wall thickness': 5.0,
                'Ostium inter-vessel distance': 0.0,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 5.0,
                'Vessel wall thickness': 5.0,
                'Vessel angle 1 degrees': 0.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': 0.0,
                'Use linear through vessel wall': True,
                'Use cross derivatives': False,
                'Refine': False,
                'Refine number of elements around': 4,
                'Refine number of elements along': 4,
                'Refine number of elements through wall': 1
            },
        }),
        'Mouse 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 12,
                'Number of elements along': 2,
                'Number of elements through wall': 1,  # not implemented for > 1
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 1.5,
                'Ostium length': 1.5,
                'Ostium wall thickness': 0.25,
                'Ostium inter-vessel distance': 0.0,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 0.5,
                'Vessel wall thickness': 0.2,
                'Vessel angle 1 degrees': 0.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': 0.0,
                'Use linear through vessel wall': True,
                'Use cross derivatives': False,
                'Refine': False,
                'Refine number of elements around': 4,
                'Refine number of elements along': 4,
                'Refine number of elements through wall': 1
            },
        }),
        'Rat 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 12,
                'Number of elements along': 2,
                'Number of elements through wall': 1,  # not implemented for > 1
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 4.0,
                'Ostium length': 3.5,
                'Ostium wall thickness': 0.5,
                'Ostium inter-vessel distance': 0.0,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 1.25,
                'Vessel wall thickness': 0.5,
                'Vessel angle 1 degrees': 0.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': 0.0,
                'Use linear through vessel wall': True,
                'Use cross derivatives': False,
                'Refine': False,
                'Refine number of elements around': 4,
                'Refine number of elements along': 4,
                'Refine number of elements through wall': 1
            },
        }),
    }

    @staticmethod
    def getName():
        return '3D Stomach 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Mouse 1',
            'Rat 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Mouse 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Mouse 1']
        elif 'Rat 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Rat 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Rat 1']
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Human 1']

        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements around esophagus': 8,
            'Number of elements around duodenum': 12,
            'Number of elements between cardia and duodenum': 6,
            'Number of elements across cardia': 1,
            'Wall thickness': 5.0,
            'Limiting ridge': False,
            'Gastro-esophagal junction': copy.deepcopy(ostiumOption),
            'Gastro-esophagal junction position along factor': 0.35,
            'Cardia derivative factor': 1.0,
            'Use linear through wall' : True,
            'Refine': False,
            'Refine number of elements surface': 4,
            'Refine number of elements through wall': 1
        }
        if 'Mouse 1' in parameterSetName:
            options['Number of elements around esophagus'] = 12
            options['Number of elements around duodenum'] = 14
            options['Number of elements between cardia and duodenum'] = 2
            options['Wall thickness'] = 0.25
            options['Gastro-esophagal junction position along factor'] = 0.53
            options['Cardia derivative factor'] = 0.3
            options['Limiting ridge'] = True
        elif 'Rat 1' in parameterSetName:
            options['Number of elements around esophagus'] = 12
            options['Number of elements around duodenum'] = 14
            options['Number of elements between cardia and duodenum'] = 2
            options['Wall thickness'] = 0.5
            options['Gastro-esophagal junction position along factor'] = 0.55
            options['Cardia derivative factor'] = 0.2
            options['Limiting ridge'] = True

        cls.updateSubScaffoldOptions(options)
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements around esophagus',
            'Number of elements around duodenum',
            'Number of elements between cardia and duodenum',
            'Number of elements across cardia',
            'Wall thickness',
            'Limiting ridge',
            'Gastro-esophagal junction',
            'Gastro-esophagal junction position along factor',
            'Cardia derivative factor',
            'Use linear through wall',
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through wall']

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [MeshType_1d_path1]
        if optionName == 'Gastro-esophagal junction':
            return [MeshType_3d_ostium1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        if optionName == 'Gastro-esophagal junction':
            return list(cls.ostiumDefaultScaffoldPackages.keys())
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
        if optionName == 'Gastro-esophagal junction':
            if not parameterSetName:
                parameterSetName = list(cls.ostiumDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.ostiumDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        if not options['Gastro-esophagal junction'].getScaffoldType() in cls.getOptionValidScaffoldTypes(
                'Gastro-esophagal junction'):
            options['Gastro-esophagal junction'] = cls.getOptionScaffoldPackage('Gastro-esophagal junction',
                                                                                 MeshType_3d_ostium1)
        if options['Number of elements around esophagus'] < 8:
            options['Number of elements around esophagus'] = 8
        if options['Number of elements around esophagus'] % 4 > 0:
            options['Number of elements around esophagus'] = options['Number of elements around esophagus'] // 4 * 4
        if options['Number of elements around duodenum'] < 12:
            options['Number of elements around duodenum'] = 12
        if options['Number of elements around duodenum'] % 2:
            options['Number of elements around duodenum'] += 1
        if options['Number of elements between cardia and duodenum'] < 2:
            options['Number of elements between cardia and duodenum'] = 2
        if options['Cardia derivative factor'] <= 0.0:
            options['Cardia derivative factor'] = 0.1
        for key in [
            'Refine number of elements surface',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1

        cls.updateSubScaffoldOptions(options)

    @classmethod
    def updateSubScaffoldOptions(cls, options):
        '''
        Update ostium sub-scaffold options which depend on parent options.
        '''
        wallThickness = options['Wall thickness']
        ostiumOptions = options['Gastro-esophagal junction']
        ostiumSettings = ostiumOptions.getScaffoldSettings()
        ostiumSettings['Ostium wall thickness'] = wallThickness

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        cls.updateSubScaffoldOptions(options)
        centralPath = options['Central path']
        elementsCountAroundEso = options['Number of elements around esophagus']
        elementsCountAroundDuod = options['Number of elements around duodenum']
        elementsAlongCardiaToDuod = options['Number of elements between cardia and duodenum']
        elementsCountThroughWall = 1
        wallThickness = options['Wall thickness']
        useCrossDerivatives = False
        useCubicHermiteThroughWall = not (options['Use linear through wall'])

        GEJPositionAlongFactor = options['Gastro-esophagal junction position along factor']
        GEJOptions = options['Gastro-esophagal junction']
        GEJSettings = GEJOptions.getScaffoldSettings()
        elementsAlongEsophagus = GEJSettings['Number of elements along']
        elementsThroughEsophagusWall = GEJSettings['Number of elements through wall']
        limitingRidge = options['Limiting ridge']
        elementsCountAcrossCardia = options['Number of elements across cardia']
        cardiaDerivativeFactor = options['Cardia derivative factor']

        elementsCountAlongTrackSurface = 20
        elementsAroundHalfEso = int(elementsCountAroundEso * 0.5)
        elementsAroundQuarterEso = int(elementsCountAroundEso * 0.25)
        elementsAroundHalfDuod = int(elementsCountAroundDuod * 0.5)
        elementsCountAlong = elementsAlongCardiaToDuod + elementsAroundHalfEso + 1
        zero = [0.0, 0.0, 0.0]

        fm = region.getFieldmodule()
        fm.beginChange()

        coordinates = findOrCreateFieldCoordinates(fm)
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        if useCubicHermiteThroughWall:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
            if useCrossDerivatives:
                nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
                nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
                nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        cache = fm.createFieldcache()
        mesh = fm.findMeshByDimension(3)

        nodeIdentifier = 1
        elementIdentifier = 1

        # Extract length of each group along stomach from central path
        arcLengthOfGroupsAlong = []
        stomachTermsAlong = [None, 'fundus of stomach', 'body of stomach', 'pyloric antrum', 'pylorus', 'duodenum']
        for i in range(len(stomachTermsAlong)):
            tmpRegion = region.createRegion()
            centralPath.generate(tmpRegion)
            cxGroup, cd1Group, cd2Group, cd3Group, cd12Group, cd13Group = \
                extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                            Node.VALUE_LABEL_D_DS2,  Node.VALUE_LABEL_D_DS3,
                                                            Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3],
                                                groupName=stomachTermsAlong[i])

            arcLength = 0.0
            for e in range(len(cxGroup) - 1):
                arcLength += interp.getCubicHermiteArcLength(cxGroup[e], cd1Group[e],
                                                             cxGroup[e + 1], cd1Group[e + 1])
            arcLengthOfGroupsAlong.append(arcLength)

            if i == 0:
                cx = cxGroup
                cd1 = cd1Group
                cd2 = cd2Group
                cd3 = cd3Group
                cd12 = cd12Group
                cd13 = cd13Group
            del tmpRegion

        sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlongTrackSurface)
        sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)
        sd3, sd13 = interp.interpolateSampleCubicHermite(cd3, cd13, se, sxi, ssf)

        # Create annotation groups for stomach sections
        allAnnotationGroups = []
        stomachCentralPathLength = sum(arcLengthOfGroupsAlong[1:])
        fundusEndPositionAlongFactor = arcLengthOfGroupsAlong[1]/stomachCentralPathLength
        arcLengthRatioForGroupsFromFundusEnd = []
        for i in range(2, len(stomachTermsAlong)):
            arcLengthRatio = (arcLengthOfGroupsAlong[i])/(stomachCentralPathLength - arcLengthOfGroupsAlong[1])
            arcLengthRatioForGroupsFromFundusEnd.append(arcLengthRatio)

        stomachGroup = AnnotationGroup(region, get_stomach_term("stomach"))
        fundusGroup = AnnotationGroup(region, get_stomach_term("fundus of stomach"))
        bodyGroup = AnnotationGroup(region, get_stomach_term("body of stomach"))
        antrumGroup = AnnotationGroup(region, get_stomach_term("pyloric antrum"))
        pylorusGroup = AnnotationGroup(region, get_stomach_term("pylorus"))
        duodenumGroup = AnnotationGroup(region, get_stomach_term("duodenum"))

        annotationGroupAlong = [[stomachGroup, fundusGroup],
                                [stomachGroup, bodyGroup],
                                [stomachGroup, antrumGroup],
                                [stomachGroup, pylorusGroup],
                                [stomachGroup, duodenumGroup]]

        # annotation fiducial points
        markerGroup = findOrCreateFieldGroup(fm, "marker")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)

        # Fundus diameter
        fundusRadius = vector.magnitude(sd2[0])
        lengthElementAlongTrackSurface = stomachCentralPathLength / elementsCountAlongTrackSurface
        elementsAlongFundus = int(fundusRadius / lengthElementAlongTrackSurface)

        d2Apex = []
        d2 = sd2[0]
        for n1 in range(elementsCountAroundDuod):
            rotAngle = n1 * 2.0 * math.pi / elementsCountAroundDuod
            rotAxis = vector.normalise(vector.crossproduct3(vector.normalise(sd2[0]), vector.normalise(sd3[0])))
            rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, rotAngle)
            d2Rot = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)]
            d2Apex.append(d2Rot)

        xEllipses = []
        d1Ellipses = []
        for n in range(elementsAlongFundus + 1, len(sx)):
            px, pd1 = createEllipsePoints(sx[n], 2 * math.pi, sd2[n], sd3[n], elementsCountAroundDuod,
                                          startRadians=0.0)
            xEllipses.append(px)
            d1Ellipses.append(pd1)

        # Find d2
        d2Raw = []
        for n1 in range(elementsCountAroundDuod):
            xAlong = []
            d2Along = []
            for n2 in range(len(xEllipses) - 1):
                v1 = xEllipses[n2][n1]
                v2 = xEllipses[n2 + 1][n1]
                d2 = findDerivativeBetweenPoints(v1, v2)
                xAlong.append(v1)
                d2Along.append(d2)
            xAlong.append(xEllipses[-1][n1])
            d2Along.append(d2)
            d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along)
            d2Raw.append(d2Smoothed)

        # Rearrange d2
        d2Ellipses = []
        for n2 in range(len(xEllipses)):
            d2Around = []
            for n1 in range(elementsCountAroundDuod):
                d2 = d2Raw[n1][n2]
                d2Around.append(d2)
            d2Ellipses.append(d2Around)

        # Merge fundus apex and body
        xAll = [[sx[0]] * elementsCountAroundDuod] + xEllipses
        d2All = [d2Apex] + d2Ellipses

        # Spread out elements
        xRaw = []
        d2Raw = []
        for n1 in range(elementsCountAroundDuod):
            xAlong = []
            d2Along = []
            for n2 in range(len(xAll)):
                xAlong.append(xAll[n2][n1])
                d2Along.append(d2All[n2][n1])
            xSampledAlong, d2SampledAlong = interp.sampleCubicHermiteCurves(xAlong, d2Along,
                                                                            elementsCountAlongTrackSurface,
                                                                            arcLengthDerivatives=True)[0:2]
            d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xSampledAlong, d2SampledAlong)
            xRaw.append(xSampledAlong)
            d2Raw.append(d2Smoothed)

        # Rearrange x and d2
        xSampledAll = []
        d1SampledAll = []
        d2SampledAll = []
        for n2 in range(elementsCountAlongTrackSurface + 1):
            xAround = []
            d1Around = []
            d2Around = []
            for n1 in range(elementsCountAroundDuod):
                x = xRaw[n1][n2]
                d2 = d2Raw[n1][n2]
                xAround.append(x)
                d2Around.append(d2)

                # Calculate d1
                if n2 > 0:
                    v1 = xRaw[n1][n2]
                    v2 = xRaw[n1 + 1 if n1 < elementsCountAroundDuod - 2 else 0][n2]
                    d1 = findDerivativeBetweenPoints(v1, v2)
                    d1Around.append(d1)
                else:
                    d1Around.append(d2Raw[int(elementsCountAroundDuod * 0.75)][0])

            if n2 > 0:
                d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)
            else:
                d1Smoothed = d1Around

            xSampledAll.append(xAround)
            d1SampledAll.append(d1Smoothed)
            d2SampledAll.append(d2Around)

        # Create tracksurface
        xTrackSurface = []
        d1TrackSurface = []
        d2TrackSurface = []
        for n2 in range(elementsCountAlongTrackSurface + 1):
            for n1 in range(elementsCountAroundDuod):
                xTrackSurface.append(xSampledAll[n2][n1])
                d1TrackSurface.append(d1SampledAll[n2][n1])
                d2TrackSurface.append(d2SampledAll[n2][n1])

        trackSurfaceStomach = TrackSurface(elementsCountAroundDuod, elementsCountAlongTrackSurface,
                                           xTrackSurface, d1TrackSurface, d2TrackSurface, loop1=True)

        # Set up gastro-esophagal junction
        GEJSettings['Number of elements around ostium'] = elementsCountAroundEso
        GEJPosition = trackSurfaceStomach.createPositionProportion(0.5, GEJPositionAlongFactor)
        xCentre, d1Centre, d2Centre = trackSurfaceStomach.evaluateCoordinates(GEJPosition, derivatives=True)
        axis1 = d1Centre

        esophagusGroup = AnnotationGroup(region, get_stomach_term("esophagus"))
        esophagusMeshGroup = esophagusGroup.getMeshGroup(mesh)
        esophagogastricJunctionGroup = AnnotationGroup(region, get_stomach_term("esophagogastric junction"))
        esophagogastricJunctionMeshGroup = esophagogastricJunctionGroup.getMeshGroup(mesh)
        stomachMeshGroup = stomachGroup.getMeshGroup(mesh)
        allAnnotationGroups += [esophagusGroup, esophagogastricJunctionGroup]

        nextNodeIdentifier, nextElementIdentifier, (o1_x, o1_d1, o1_d2, o1_d3, o1_NodeId, o1_Positions) = \
            generateOstiumMesh(region, GEJSettings, trackSurfaceStomach, GEJPosition, axis1,
                               nodeIdentifier, elementIdentifier, vesselMeshGroups=[[stomachMeshGroup, esophagusMeshGroup]],
                               ostiumMeshGroups=[stomachMeshGroup, esophagogastricJunctionMeshGroup])

        stomachStartNode = nextNodeIdentifier
        stomachStartElement = nextElementIdentifier
        nodeIdentifier = nextNodeIdentifier
        elementIdentifier = nextElementIdentifier

        fundusEndPosition = trackSurfaceStomach.createPositionProportion(0.0, fundusEndPositionAlongFactor)
        xFundusEnd, d1FundusEnd, d2FundusEnd = trackSurfaceStomach.evaluateCoordinates(fundusEndPosition, derivatives=True)
        elementsAlongFundus = elementsAroundQuarterEso + (0 if limitingRidge else 1)

        # From fundus end to duodenum
        elementsAlongFundusEndToDuod = elementsCountAlong - elementsAlongFundus
        xAlongDownFundus = []
        d2AlongDownFundus = []
        xAlongDownFundus.append(xFundusEnd)
        d2AlongDownFundus.append(d2FundusEnd)
        startIdx = fundusEndPositionAlongFactor*elementsCountAlongTrackSurface
        startIdx = math.ceil(startIdx) + (1 if startIdx - math.ceil(startIdx) == 0 else 0)
        for n2 in range(startIdx, len(xSampledAll)):
            xAlongDownFundus.append(xSampledAll[n2][0])
            d2AlongDownFundus.append(d2SampledAll[n2][0])

        # From esophagus to fundus end
        elementsAlongGCFromEsoToFundusEnd = elementsAroundHalfDuod - 2 + elementsAlongFundus
        xAlongUpFundus = []
        d2AlongUpFundus = []
        xAlongUpFundus.append(o1_x[1][0])
        d2AlongUpFundus.append([cardiaDerivativeFactor * c for c in o1_d2[1][0]])

        # From esophagus to fundus apex
        esoStartProportion2 = trackSurfaceStomach.getProportion(o1_Positions[0])[1]
        elementsAlongUpstreamOfEso = int(elementsCountAlongTrackSurface * esoStartProportion2)
        arcLengthEsoApex = 0.0
        for n2 in range(elementsAlongUpstreamOfEso):
            nAlong = elementsAlongUpstreamOfEso - n2
            v1 = xSampledAll[nAlong][elementsAroundHalfDuod]
            v2 = xSampledAll[nAlong - 1][elementsAroundHalfDuod]
            d = [v2[c] - v1[c] for c in range(3)]
            arcLengthAround = interp.computeCubicHermiteArcLength(v1, d, v2, d, True)
            arcLengthEsoApex += arcLengthAround
            d2 = [c * arcLengthAround for c in vector.normalise(d)]
            xAlongUpFundus.append(v1)
            d2AlongUpFundus.append(d2)

        # From fundus apex to end
        for n2 in range(startIdx + (1 if limitingRidge else 0)):
            xAlongUpFundus.append(xSampledAll[n2][0])
            d2AlongUpFundus.append(d2SampledAll[n2][0])

        if limitingRidge:
            # Sample from esophagus to fundus end
            xAlongGCEsoToFundusEnd, d2AlongGCEsoToFundusEnd = \
                interp.sampleCubicHermiteCurvesSmooth(xAlongUpFundus, d2AlongUpFundus,
                                                      elementsAlongGCFromEsoToFundusEnd,
                                                      derivativeMagnitudeStart=vector.magnitude(d2AlongUpFundus[0]))[0:2]

            # Sample from limiting ridge to duodenum
            xAlongDownFundus[0] = xAlongGCEsoToFundusEnd[-1]
            d2AlongDownFundus[0] = d2AlongGCEsoToFundusEnd[-1]

            if xAlongDownFundus[0] == xAlongDownFundus[1]:
                del xAlongDownFundus[1], d2AlongDownFundus[1]

            xAlongGCFundusEndToDuod, d2AlongGCFundusEndToDuod = \
                interp.sampleCubicHermiteCurves(xAlongDownFundus, d2AlongDownFundus, elementsAlongFundusEndToDuod,
                                                addLengthStart=0.5 * vector.magnitude(d2AlongDownFundus[0]),
                                                lengthFractionStart=0.5, arcLengthDerivatives=True)[0:2]

        else:
            # Sample from end of fundus to duodenum
            xAlongGCFundusEndToDuod, d2AlongGCFundusEndToDuod = \
                interp.sampleCubicHermiteCurves(xAlongDownFundus, d2AlongDownFundus, elementsAlongFundusEndToDuod,
                                                arcLengthDerivatives=True)[0:2]
            # Sample from ostium to end of fundus
            # Provide first element size for transition
            xAlongUpFundus[-1] = xAlongGCFundusEndToDuod[0]
            d2AlongUpFundus[-1] = d2AlongGCFundusEndToDuod[0]
            xAlongGCEsoToFundusEnd, d2AlongGCEsoToFundusEnd = \
                interp.sampleCubicHermiteCurvesSmooth(xAlongUpFundus, d2AlongUpFundus,
                                                      elementsAlongGCFromEsoToFundusEnd,
                                                      derivativeMagnitudeStart= vector.magnitude(d2AlongUpFundus[0]),
                                                      derivativeMagnitudeEnd = vector.magnitude(d2AlongUpFundus[-1]))[0:2]

        xEsoToDuodGC = xAlongGCEsoToFundusEnd[:-1] + xAlongGCFundusEndToDuod
        d2EsoToDuodGC = d2AlongGCEsoToFundusEnd[:-1] + d2AlongGCFundusEndToDuod

        # Use element spacing on greater curvature to decide where annotation groups should sit
        arcLengthGCFromFundusEnd = 0.0
        arcLengthsGCFromFundusEnd = []
        for n2 in range(len(xAlongGCFundusEndToDuod) - 1):
            arcLengthGCFromFundusEnd += interp.getCubicHermiteArcLength(xAlongGCFundusEndToDuod[n2],
                                                                        d2AlongGCFundusEndToDuod[n2],
                                                                        xAlongGCFundusEndToDuod[n2 + 1],
                                                                        d2AlongGCFundusEndToDuod[n2 + 1])
            arcLengthsGCFromFundusEnd.append(arcLengthGCFromFundusEnd)
        arcLengthGCFundusEndToDuod = arcLengthGCFromFundusEnd

        # Move along elements on GC and assign to annotation group
        elementsCountAlongGroups = []
        elementsCountAlongGroups.append(elementsAlongFundus)
        groupLength = 0.0
        e = 0
        elementsCount = 1
        length = arcLengthsGCFromFundusEnd[e]
        for i in range(len(arcLengthRatioForGroupsFromFundusEnd)):
            groupLength += arcLengthRatioForGroupsFromFundusEnd[i] * arcLengthGCFundusEndToDuod
            if e == len(arcLengthsGCFromFundusEnd) - 2:
                elementsCount += 1
                elementsCountAlongGroups.append(elementsCount)
            else:
                while length < groupLength:
                    elementsCount += 1
                    e += 1
                    length = arcLengthsGCFromFundusEnd[e]

                # check which end is grouplength closer to
                distToUpperEnd = abs(length - groupLength)
                distToLowerEnd = abs(groupLength - arcLengthsGCFromFundusEnd[e - 1])
                if distToLowerEnd < distToUpperEnd:
                    elementsCount -= 1
                    elementsCountAlongGroups.append(elementsCount)
                    e -= 1
                    length = arcLengthsGCFromFundusEnd[e]
                else:
                    elementsCountAlongGroups.append(elementsCount)
            elementsCount = 0

        annotationGroupsAlong = []
        for i in range(len(elementsCountAlongGroups)):
            elementsCount = elementsCountAlongGroups[i]
            for n in range(elementsCount):
                annotationGroupsAlong.append(annotationGroupAlong[i])

        arcLength = 0.0
        for e in range(len(xEsoToDuodGC) - 1):
            arcLength += interp.getCubicHermiteArcLength(xEsoToDuodGC[e], d2EsoToDuodGC[e],
                                                         xEsoToDuodGC[e + 1], d2EsoToDuodGC[e + 1])
            if arcLength > arcLengthEsoApex:
                nodesCountFromEsoToApex = e + 2
                break

        ptsOnTrackSurfaceGC = []
        xAroundEllipse = []
        d1AroundEllipse = []
        for n2 in range(elementsCountAlongTrackSurface + 1):
            ptsOnTrackSurfaceGC.append(xSampledAll[n2][0])

        # Ring upstream of 6 point junction
        # First half
        n2 = elementsAroundHalfEso - 1
        ostiumIdx = n2
        GCIdx = elementsAroundHalfDuod - 1 + n2
        GCPosition, d1GC = findClosestPositionAndDerivativeOnTrackSurface(xEsoToDuodGC[GCIdx], ptsOnTrackSurfaceGC,
                                                                          trackSurfaceStomach, 0.0,
                                                                          elementsCountAlongTrackSurface)
        GCProportion1, GCProportion2 = trackSurfaceStomach.getProportion(GCPosition)
        endPosition = o1_Positions[ostiumIdx]
        rotFrame = matrix.getRotationMatrixFromAxisAngle(vector.normalise(o1_d1[1][ostiumIdx]), math.pi)
        d2 = o1_d2[1][ostiumIdx]
        d1EndOstium = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)]
        endProportion1, endProportion2 = trackSurfaceStomach.getProportion(endPosition)

        xFirstHalf = \
            getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, 0.0, GCProportion2, endProportion1,
                                                   endProportion2, elementsAroundHalfDuod + 1,
                                                   startDerivative=d1GC, endDerivative=d1EndOstium,
                                                   endDerivativeMagnitude=cardiaDerivativeFactor * vector.magnitude(d2))[0]
        # Second half
        ostiumIdx2 = -n2
        startPosition = o1_Positions[ostiumIdx2]
        d1StartOstium = o1_d2[1][ostiumIdx2]
        startProportion1, startProportion2 = trackSurfaceStomach.getProportion(startPosition)
        xSecondHalf = \
            getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, startProportion1, startProportion2, 1.0,
                                                   GCProportion2, elementsAroundHalfDuod + 1,
                                                   startDerivative=d1StartOstium, endDerivative=d1GC,
                                                   startDerivativeMagnitude=cardiaDerivativeFactor *
                                                                            vector.magnitude(d1StartOstium))[0]

        xAroundBefore6Pt = xFirstHalf[:-1] + xSecondHalf[1:-1]

        # Rings downstream of 6 pt junction
        for idx in range(-(elementsCountAlong - elementsAroundHalfEso - 1), 0):
            # Search for point on central path and use that to make ellipse
            xStart = xEsoToDuodGC[idx]
            startPosition, d1Start = findClosestPositionAndDerivativeOnTrackSurface(xEsoToDuodGC[idx],
                                                                                    ptsOnTrackSurfaceGC,
                                                                                    trackSurfaceStomach, 0.0,
                                                                                    elementsCountAlongTrackSurface)
            startProportion2 = trackSurfaceStomach.getProportion(startPosition)[1]
            if 0.0 < startProportion2 < 1.0:
                closestIdxOnCentralPath = interp.getNearestPointIndex(sx, xStart)
                if 0 < closestIdxOnCentralPath < len(sx) - 1:
                    # Check if xStart is closer to upstream or downstream of closestIdx
                    xOnGCPrevElem = [sx[closestIdxOnCentralPath - 1][c] + sd2[closestIdxOnCentralPath - 1][c] for c in range(3)]
                    distBetweenXOnGCPrevElem = vector.magnitude([xStart[c] - xOnGCPrevElem[c] for c in range(3)])
                    xOnGCNextElem = [sx[closestIdxOnCentralPath + 1][c] + sd2[closestIdxOnCentralPath + 1][c] for c in range(3)]
                    distBetweenXOnGCNextElem = vector.magnitude([xStart[c] - xOnGCNextElem[c] for c in range(3)])
                    eiLowerLimit = closestIdxOnCentralPath - (1 if distBetweenXOnGCNextElem > distBetweenXOnGCPrevElem else 0)
                elif closestIdxOnCentralPath == len(sx) - 1:
                    eiLowerLimit = closestIdxOnCentralPath - 1
                elif closestIdxOnCentralPath == 0:
                    eiLowerLimit = closestIdxOnCentralPath

                xiLowerLimit = 0.0
                xiUpperLimit = 1.0
                xOnGCLowerLimit = [sx[eiLowerLimit][c] + sd2[eiLowerLimit][c] for c in range(3)]
                xOnGCUpperLimit = [sx[eiLowerLimit + 1][c] + sd2[eiLowerLimit + 1][c] for c in range(3)]
                distBetweenXAndXStartPrev = vector.magnitude([xStart[c] - xOnGCLowerLimit[c] for c in range(3)])

                # Search for point on central path which is orthogonal to xStart
                tol = 1e-8
                xLowerLimit = sx[eiLowerLimit]
                d1LowerLimit = sd1[eiLowerLimit]
                d2LowerLimit = sd2[eiLowerLimit]
                d12LowerLimit = sd12[eiLowerLimit]
                xUpperLimit = sx[eiLowerLimit + 1]
                d1UpperLimit = sd1[eiLowerLimit + 1]
                d2UpperLimit = sd2[eiLowerLimit + 1]
                d12UpperLimit = sd12[eiLowerLimit + 1]

                for iter in range(100):
                    xiGuess = 0.5 * (xiLowerLimit + xiUpperLimit)
                    x = interp.interpolateCubicHermite(xLowerLimit, d1LowerLimit, xUpperLimit, d1UpperLimit, xiGuess)
                    d2 = interp.interpolateCubicHermite(d2LowerLimit, d12LowerLimit, d2UpperLimit, d12UpperLimit, xiGuess)
                    xGuess = [x[c] + d2[c] for c in range(3)]
                    distBetweenXAndXStart = vector.magnitude([xStart[c] - xGuess[c] for c in range(3)])
                    distBetweenXOnGCLowerLimitAndXStart = vector.magnitude([xStart[c] - xOnGCLowerLimit[c] for c in range(3)])
                    distBetweenXOnGCUpperLimitAndXStart = vector.magnitude([xStart[c] - xOnGCUpperLimit[c] for c in range(3)])

                    if abs(distBetweenXAndXStart - distBetweenXAndXStartPrev) < tol:
                        xProjection = x
                        xiProjection = xiGuess
                        break
                    elif distBetweenXOnGCLowerLimitAndXStart < distBetweenXOnGCUpperLimitAndXStart:
                        xiUpperLimit = xiGuess
                        xOnGCUpperLimit = xGuess
                        distBetweenXAndXStartPrev = distBetweenXAndXStart
                    else:
                        xiLowerLimit = xiGuess
                        xOnGCLowerLimit = xGuess
                        distBetweenXAndXStartPrev = distBetweenXAndXStart

                if iter > 98:
                    print('Search for projection on central path - Max iters reached:', iter)

                axis1 = interp.interpolateCubicHermite(sd2[eiLowerLimit], sd12[eiLowerLimit],
                                                       sd2[eiLowerLimit + 1], sd12[eiLowerLimit + 1], xiProjection)
                axis2 = interp.interpolateCubicHermite(sd3[eiLowerLimit], sd13[eiLowerLimit],
                                                       sd3[eiLowerLimit + 1], sd13[eiLowerLimit + 1], xiProjection)
                xAround, d1Around = createEllipsePoints(xProjection, 2 * math.pi, axis1, axis2, elementsCountAroundDuod,
                                                        startRadians=0.0)

            elif startProportion2 == 0:
                xAround, d1Around = createEllipsePoints(sx[0], 2 * math.pi, sd2[0], sd3[0], elementsCountAroundDuod,
                                                        startRadians=0.0)
            elif startProportion2 == 1.0:
                xAround, d1Around = createEllipsePoints(sx[-1], 2 * math.pi, sd2[-1], sd3[-1], elementsCountAroundDuod,
                                                        startRadians=0.0)

            d1Around = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around,
                                                                magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
            xAroundEllipse.append(xAround)
            d1AroundEllipse.append(d1Around)

            # Average adjacent ring with first downstream ring that is not adjacent to esophagus
            if idx == -(elementsCountAlong - elementsAroundHalfEso - 1):
                xAve = []
                xAve.append(xEsoToDuodGC[idx - 1])

                for n in range(1, elementsCountAroundDuod):
                    startPosition = trackSurfaceStomach.findNearestPosition(xAround[n])
                    startProportion1, startProportion2 = trackSurfaceStomach.getProportion(startPosition)
                    if n == elementsAroundHalfDuod:
                        endPosition = o1_Positions[elementsAroundHalfEso]
                    else:
                        endPosition = trackSurfaceStomach.findNearestPosition(xAroundBefore6Pt[n + (0 if n < elementsAroundHalfDuod else 1)])
                    endProportion1, endProportion2 = trackSurfaceStomach.getProportion(endPosition)
                    xSampled = getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, startProportion1,
                                                                      startProportion2, endProportion1, endProportion2, 2)[0]
                    xAve.append(xSampled[1])

                # Find 6 pt junction
                p1x_6pt = o1_x[1][elementsAroundHalfEso]
                d = o1_d2[1][elementsAroundHalfEso]
                rotFrame = matrix.getRotationMatrixFromAxisAngle(o1_d1[1][elementsAroundHalfEso], math.pi)
                p1d = [rotFrame[j][0] * d[0] + rotFrame[j][1] * d[1] + rotFrame[j][2] * d[2] for j in range(3)]
                p1d_6pt = [cardiaDerivativeFactor * c for c in p1d]

                p2x_6pt = xAve[int(len(xAve)*0.5) + 1]
                p2d_6pt = findDerivativeBetweenPoints(p2x_6pt, xAve[int(len(xAve) * 0.5) + 2])

                p3x_6pt = xAve[int(len(xAve) * 0.5) - 1]
                p3d_6pt = findDerivativeBetweenPoints(p3x_6pt, xAve[int(len(xAve) * 0.5) - 2])

                x6pt = get_bifurcation_triple_point(p1x_6pt, p1d_6pt, p2x_6pt, p2d_6pt, p3x_6pt, p3d_6pt)[0]

        # Gradually vary derivative magnitude at annulus to create smooth transition of nodes around the cardia
        # between triple point to 6pt junction
        distBetween6ptJunctionOstium = vector.magnitude([o1_x[1][elementsAroundHalfEso][c] - x6pt[c] for c in range(3)])
        distAnnulusAtQuarterEso = cardiaDerivativeFactor * vector.magnitude(o1_d2[1][elementsAroundQuarterEso])

        n = 0
        xAlongAround = []
        d1AlongAround = []
        for n2 in range(elementsAroundQuarterEso + 1, elementsAroundHalfEso):
            xi = n/elementsAroundQuarterEso
            derivativeMagnitude = xi * distBetween6ptJunctionOstium + (1-xi) * distAnnulusAtQuarterEso
            ostiumIdx = n2
            GCIdx = elementsAroundHalfDuod - 1 + n2
            GCPosition, d1GC = findClosestPositionAndDerivativeOnTrackSurface(xEsoToDuodGC[GCIdx], ptsOnTrackSurfaceGC,
                                                                              trackSurfaceStomach, 0.0,
                                                                              elementsCountAlongTrackSurface)
            GCProportion1, GCProportion2 = trackSurfaceStomach.getProportion(GCPosition)
            endPosition = o1_Positions[ostiumIdx]
            rotFrame = matrix.getRotationMatrixFromAxisAngle(vector.normalise(o1_d1[1][ostiumIdx]), math.pi)
            d2 = o1_d2[1][ostiumIdx]
            d1EndOstium = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)]
            endProportion1, endProportion2 = trackSurfaceStomach.getProportion(endPosition)

            xFirstHalf, d1FirstHalf = \
                getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, 0.0, GCProportion2, endProportion1,
                                                       endProportion2, elementsAroundHalfDuod + 1,
                                                       startDerivative=d1GC, endDerivative=d1EndOstium,
                                                       endDerivativeMagnitude=derivativeMagnitude)
            # Second half
            ostiumIdx2 = -n2
            startPosition = o1_Positions[ostiumIdx2]
            d1StartOstium = o1_d2[1][ostiumIdx2]
            startProportion1, startProportion2 = trackSurfaceStomach.getProportion(startPosition)
            xSecondHalf, d1SecondHalf = \
                getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, startProportion1, startProportion2, 1.0,
                                                       GCProportion2, elementsAroundHalfDuod + 1,
                                                       startDerivative=d1StartOstium, endDerivative=d1GC,
                                                       startDerivativeMagnitude=derivativeMagnitude)

            xAround = xFirstHalf[:-1] + xSecondHalf[1:-1]
            d1Around = d1FirstHalf[:-1] + d1SecondHalf[1:-1]

            n += 1
            xAlongAround.append(xAround)
            d1AlongAround.append(d1Around)

        # Resample ring with 6pt junction to improve smoothness
        idx = -(elementsCountAlong - elementsAroundHalfEso - 1)
        xAve = []
        dAve = []
        xAve.append(xEsoToDuodGC[idx - 1])

        for n1 in range(1, elementsCountAroundDuod + 1):
            startPosition = trackSurfaceStomach.findNearestPosition(xAlongAround[-1][n1])
            startProportion1, startProportion2 = trackSurfaceStomach.getProportion(startPosition)
            endPosition = trackSurfaceStomach.findNearestPosition(xAroundEllipse[0][n1 + (0 if n1 < elementsAroundHalfDuod + 1 else -1)])
            endProportion1, endProportion2 = trackSurfaceStomach.getProportion(endPosition)
            xSampled = getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, startProportion1,
                                                              startProportion2, endProportion1, endProportion2, 2)[0]
            xAve.append(xSampled[1])

        xAve[int(len(xAve) * 0.5)] = x6pt
        del xAve[int(len(xAve) * 0.5) + 1]

        for n1 in range(len(xAve)):
            v1 = xAve[n1]
            v2 = xAve[(n1 + 1) % len(xAve)]
            d1 = findDerivativeBetweenPoints(v1, v2)
            dAve.append(d1)
        dAve = interp.smoothCubicHermiteDerivativesLoop(xAve, dAve)
        xAlongAround.append(xAve)
        d1AlongAround.append(dAve)

        xAlongAround += xAroundEllipse
        d1AlongAround += d1AroundEllipse

        # Sample 2 loops next to annulus from point on GC to point on first ring on xAlongAround
        ptsOnTrackSurfaceEsoToFundus = []
        for n2 in range(elementsCountAlongTrackSurface + 1):
            ptsOnTrackSurfaceEsoToFundus.append(xSampledAll[n2][elementsAroundHalfDuod])

        xLoopsRight = []
        xLoopsLeft = []
        for nLoop in range(1, elementsAroundHalfDuod - 1):
            GCIdx = nLoop + 1
            if GCIdx < nodesCountFromEsoToApex:
                ptsOnTrackSurface = ptsOnTrackSurfaceEsoToFundus
                proportion1 = 0.5
            else:
                ptsOnTrackSurface = ptsOnTrackSurfaceGC
                proportion1 = 0.0
            d2GC = findClosestPositionAndDerivativeOnTrackSurface(xEsoToDuodGC[GCIdx],
                                                                  ptsOnTrackSurface,
                                                                  trackSurfaceStomach, proportion1,
                                                                  elementsCountAlongTrackSurface)[1]
            if GCIdx < nodesCountFromEsoToApex:
                rotFrame = matrix.getRotationMatrixFromAxisAngle(vector.normalise(d2EsoToDuodGC[GCIdx]), math.pi)
                d2GCRot = [rotFrame[j][0] * d2GC[0] + rotFrame[j][1] * d2GC[1] + rotFrame[j][2] * d2GC[2] for j in range(3)]
                d2GC = d2GCRot

            for nSide in range(2):
                if nSide == 0:
                    xEnd = xAlongAround[0][elementsAroundHalfDuod - nLoop]
                    d2End = [xAlongAround[1][elementsAroundHalfDuod - nLoop][c] -
                             xAlongAround[0][elementsAroundHalfDuod - nLoop][c] for c in range(3)]
                else:
                    rotFrame = matrix.getRotationMatrixFromAxisAngle(vector.normalise(d2EsoToDuodGC[GCIdx]), math.pi)
                    d2GCRot = [rotFrame[j][0] * d2GC[0] + rotFrame[j][1] * d2GC[1] + rotFrame[j][2] * d2GC[2] for j in range(3)]
                    d2GC = d2GCRot

                    xEnd = xAlongAround[0][elementsAroundHalfDuod + 1 + nLoop]
                    d2End = [xAlongAround[1][elementsAroundHalfDuod + nLoop + 1][c] -
                             xAlongAround[0][elementsAroundHalfDuod + (1 if elementsCountAroundEso > 8 else 2) + nLoop][c] for c in range(3)]

                nx = [xEsoToDuodGC[GCIdx], xEnd]
                nd2 = [d2GC, d2End]
                x, d2 = interp.sampleCubicHermiteCurves(nx, nd2, elementsAroundQuarterEso + 2, arcLengthDerivatives=True)[0:2]

                # Find closest sampled points onto track surface
                xProjectedPoints = []
                d2ProjectedPoints = []
                xProjectedPoints.append(xEsoToDuodGC[GCIdx])
                d2ProjectedPoints.append(d2GC)
                for n2 in range(1, len(x)):
                    projectedPosition = trackSurfaceStomach.findNearestPosition(x[n2])
                    xProjected = trackSurfaceStomach.evaluateCoordinates(projectedPosition)
                    xProjectedPoints.append(xProjected)

                for n2 in range(1, len(xProjectedPoints) - 1):
                    d2 = findDerivativeBetweenPoints(xProjectedPoints[n2], xProjectedPoints[n2 + 1])
                    d2ProjectedPoints.append(d2)
                d2ProjectedPoints.append(d2End)

                # Sample points again
                xLoop = interp.sampleCubicHermiteCurves(xProjectedPoints, d2ProjectedPoints,
                                                        elementsAroundQuarterEso + 2,
                                                        addLengthEnd=0.5 * vector.magnitude(d2ProjectedPoints[-1]),
                                                        lengthFractionEnd=0.5, arcLengthDerivatives=True)[0]
                (xLoopsRight if nSide == 0 else xLoopsLeft).append(xLoop)

        # Find triple point
        xTriplePts = [[None], [None]]  # Right, left
        d1TriplePts = [[None], [None]]
        d2TriplePts = [[None], [None]]
        d3TriplePtsNorm = [[None], [None]]

        for nSide in range(2):
            ostiumIdx = elementsAroundQuarterEso if nSide == 0 else -elementsAroundQuarterEso
            p1x = o1_x[1][ostiumIdx]
            d = o1_d2[1][ostiumIdx]
            rotFrame = matrix.getRotationMatrixFromAxisAngle(o1_d1[1][ostiumIdx], math.pi)
            p1d = [rotFrame[j][0] * d[0] + rotFrame[j][1] * d[1] + rotFrame[j][2] * d[2] for j in range(3)]
            p1d = [cardiaDerivativeFactor * c for c in p1d]

            xLoops = xLoopsRight if nSide == 0 else xLoopsLeft
            p2x = xLoops[0][elementsAroundQuarterEso + 1]  # downstream bifurcation
            p2d = findDerivativeBetweenPoints(xLoops[0][elementsAroundQuarterEso + 1],
                                              xLoops[1][elementsAroundQuarterEso + 1])

            p3x = xLoops[0][elementsAroundQuarterEso]
            p3d = findDerivativeBetweenPoints(xLoops[0][elementsAroundQuarterEso],
                                              xLoops[1][elementsAroundQuarterEso])

            xTriplePts[nSide], d1TriplePts[nSide], d2TriplePts[nSide] = get_bifurcation_triple_point(p1x, p1d,
                                                                                                     p2x, p2d,
                                                                                                     p3x, p3d)
            d3TriplePtsNorm[nSide] = vector.normalise(
                vector.crossproduct3(vector.normalise(d1TriplePts[nSide]),
                                     vector.normalise(d2TriplePts[nSide])))

            # Make sure triple point is on track surface
            triplePointPosition = trackSurfaceStomach.findNearestPosition(xTriplePts[nSide])
            xTriplePts[nSide] = trackSurfaceStomach.evaluateCoordinates(triplePointPosition)

        # Sample points from GC to bottom of loops to create nodes running on row 2
        xBifurcationRings = []
        d1BifurcationRings = []
        xUp = []
        d1Up = []
        for n2 in range(elementsAroundQuarterEso):
            xAroundRight = []
            d1AroundRight = []
            xAroundLeft = []
            d1AroundLeft = []
            loopIdx = n2 + 2
            ostiumIdx = loopIdx + (0 if n2 < elementsAroundQuarterEso - 1 else -1)
            GCIdx = elementsAlongGCFromEsoToFundusEnd - elementsAroundQuarterEso + n2 + (2 if limitingRidge else 1)
            d1GC = findClosestPositionAndDerivativeOnTrackSurface(xEsoToDuodGC[GCIdx], ptsOnTrackSurfaceGC,
                                                                  trackSurfaceStomach, 0.0,
                                                                  elementsCountAlongTrackSurface)[1]
            for nSide in range(2):
                if nSide == 0: # Right side
                    xAroundRight.append(xEsoToDuodGC[GCIdx])
                    d1AroundRight.append(d1GC)
                    xOnLastLoopRight = xLoopsRight[-1][loopIdx]
                    d1OnLastLoopRight = findDerivativeBetweenPoints(xLoopsRight[-1][loopIdx], xLoopsRight[-2][loopIdx])

                    nx = [xEsoToDuodGC[GCIdx], xOnLastLoopRight]
                    nd1 = [d1GC, d1OnLastLoopRight]
                    x, d1 = interp.sampleCubicHermiteCurves(nx, nd1, 2, arcLengthDerivatives=True)[0:2]

                    # Find closest sampled points onto track surface
                    projectedPosition = trackSurfaceStomach.findNearestPosition(x[1])
                    x[1] = trackSurfaceStomach.evaluateCoordinates(projectedPosition)
                    d1[1] = findDerivativeBetweenPoints(x[1], x[2])

                    # Sample points again
                    x, d1 = interp.sampleCubicHermiteCurves(x, d1, 2, arcLengthDerivatives=True)[0:2]
                    xAroundRight.append(x[1])
                    d1AroundRight.append(d1[1])

                    for n in range(len(xLoopsRight) - 1):
                        xAroundRight.append(xLoopsRight[-(1 + n)][loopIdx])
                        d1AroundRight.append(findDerivativeBetweenPoints(xLoopsRight[-(1 + n)][loopIdx], xLoopsRight[-(1 + n + 1)][loopIdx]))

                    if loopIdx < elementsAroundQuarterEso:  # additional elements upstream of triple point
                        xLoop = xLoopsRight[0][loopIdx]
                        xLoopPosition = trackSurfaceStomach.findNearestPosition(xLoop)
                        xLoopProportion1, xLoopProportion2 = trackSurfaceStomach.getProportion(xLoopPosition)
                        xOstium = o1_x[1][ostiumIdx]
                        ostiumPosition = o1_Positions[ostiumIdx]
                        ostiumProportion1, ostiumProportion2 = trackSurfaceStomach.getProportion(ostiumPosition)
                        d = findDerivativeBetweenPoints(xLoop, xOstium)
                        endDerivativeMag = vector.magnitude(o1_d2[1][ostiumIdx]) * cardiaDerivativeFactor
                        xSampled, dSampled = \
                            getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, xLoopProportion1,
                                                                   xLoopProportion2, ostiumProportion1,
                                                                   ostiumProportion2, 2,
                                                                   endDerivativeMagnitude=endDerivativeMag)[0:2]
                        xAroundRight += xSampled[:2]
                        d1AroundRight += dSampled[:2]

                    else: # connected to triple point
                        xAroundRight += [xLoopsRight[0][loopIdx]] + [xTriplePts[0]]
                        d1AroundRight += [findDerivativeBetweenPoints(xLoopsRight[0][loopIdx], xTriplePts[0])] + \
                                         [d1TriplePts[0]]

                else: # left side
                    if loopIdx < elementsAroundQuarterEso:  # additional elements upstream of triple point
                        xLoop = xLoopsLeft[0][loopIdx]
                        xLoopPosition = trackSurfaceStomach.findNearestPosition(xLoop)
                        xLoopProportion1, xLoopProportion2 = trackSurfaceStomach.getProportion(xLoopPosition)
                        xOstium = o1_x[1][-ostiumIdx]
                        ostiumPosition = o1_Positions[-ostiumIdx]
                        ostiumProportion1, ostiumProportion2 = trackSurfaceStomach.getProportion(ostiumPosition)
                        d = findDerivativeBetweenPoints(xOstium, xLoop)
                        startDerivativeMag = vector.magnitude(o1_d2[1][-ostiumIdx]) * cardiaDerivativeFactor
                        xSampled, dSampled = \
                            getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, ostiumProportion1,
                                                                   ostiumProportion2, xLoopProportion1,
                                                                   xLoopProportion2, 2,
                                                                   startDerivativeMagnitude=startDerivativeMag)[0:2]
                        xAroundLeft.append(xSampled[1])
                        d1AroundLeft.append(dSampled[1])
                    else:
                        xAroundLeft.append(xTriplePts[1])
                        d1AroundLeft.append(d1TriplePts[1])

                    for n in range(len(xLoopsLeft) - 1):
                        xAroundLeft.append(xLoopsLeft[n][loopIdx])
                        d1AroundLeft.append(findDerivativeBetweenPoints(xLoopsLeft[n][loopIdx], xLoopsLeft[n + 1][loopIdx]))

                    xOnLastLoopLeft = xLoopsLeft[-1][loopIdx]
                    d1OnLastLoopLeft = findDerivativeBetweenPoints(xLoopsLeft[-2][loopIdx], xLoopsLeft[-1][loopIdx])

                    nx = [xOnLastLoopLeft, xEsoToDuodGC[GCIdx]]
                    nd1 = [d1OnLastLoopLeft, d1GC]
                    x, d1 = interp.sampleCubicHermiteCurves(nx, nd1, 2, arcLengthDerivatives=True)[0:2]

                    # Find closest sampled points onto track surface
                    projectedPosition = trackSurfaceStomach.findNearestPosition(x[1])
                    x[1] = trackSurfaceStomach.evaluateCoordinates(projectedPosition)
                    d1[1] = findDerivativeBetweenPoints(x[1], x[2])

                    # Sample points again
                    x, d1 = interp.sampleCubicHermiteCurves(x, d1, 2, arcLengthDerivatives=True)[0:2]
                    xAroundLeft += [xOnLastLoopLeft] + [x[1]] + [xEsoToDuodGC[GCIdx]]
                    d1AroundLeft += [findDerivativeBetweenPoints(xOnLastLoopLeft, x[1])] + [d1[1]] + [d1GC]

            xAround = xAroundRight + xAroundLeft[:-1]
            d1Around = d1AroundRight + d1AroundLeft[:-1]
            xUp.append(xAround)
            d1Up.append(d1Around)

            if loopIdx >= elementsAroundQuarterEso:
                xBifurcationRings.append(xAround)
                d1BifurcationRings.append(d1Around)

        # Row 2
        xRow2Right = []
        d1Row2Right = []
        xRow2Left = []
        d1Row2Left = []

        for nSide in range(2):
            loopIdx = 1
            ostiumIdx = 1
            if nSide == 0:
                xRow2Right.append(xUp[0][1])
                d1Row2Right.append(findDerivativeBetweenPoints(xUp[0][1], xLoopsRight[-1][1]))
                # Append rows upwards in loops
                for n in range(len(xLoopsRight) - 1):
                    xRow2Right.append(xLoopsRight[-(1 + n)][1])
                    d1Row2Right.append(
                        findDerivativeBetweenPoints(xLoopsRight[-(1 + n)][1], xLoopsRight[-(1 + n + 1)][1]))

                xLoop = xLoopsRight[0][loopIdx]
                xLoopPosition = trackSurfaceStomach.findNearestPosition(xLoop)
                xLoopProportion1, xLoopProportion2 = trackSurfaceStomach.getProportion(xLoopPosition)
                xOstium = o1_x[1][ostiumIdx]
                ostiumPosition = o1_Positions[ostiumIdx]
                ostiumProportion1, ostiumProportion2 = trackSurfaceStomach.getProportion(ostiumPosition)
                d = findDerivativeBetweenPoints(xLoop, xOstium)
                endDerivativeMag = vector.magnitude(o1_d2[1][ostiumIdx]) * cardiaDerivativeFactor
                xSampled, dSampled = \
                    getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, xLoopProportion1,
                                                           xLoopProportion2, ostiumProportion1,
                                                           ostiumProportion2, 2,
                                                           endDerivativeMagnitude=endDerivativeMag)[0:2]
                xRow2Right += xSampled[0:2]
                d1Row2Right += [findDerivativeBetweenPoints(xSampled[0], xSampled[1])] + [dSampled[1]]

            else:
                xLoop = xLoopsLeft[0][loopIdx]
                xLoopPosition = trackSurfaceStomach.findNearestPosition(xLoop)
                xLoopProportion1, xLoopProportion2 = trackSurfaceStomach.getProportion(xLoopPosition)
                xOstium = o1_x[1][-ostiumIdx]
                ostiumPosition = o1_Positions[-ostiumIdx]
                ostiumProportion1, ostiumProportion2 = trackSurfaceStomach.getProportion(ostiumPosition)
                d = findDerivativeBetweenPoints(xOstium, xLoop)
                startDerivativeMag = vector.magnitude(o1_d2[1][-ostiumIdx]) * cardiaDerivativeFactor

                xSampled, dSampled = \
                    getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, ostiumProportion1, ostiumProportion2,
                                                           xLoopProportion1, xLoopProportion2, 2,
                                                           startDerivativeMagnitude=startDerivativeMag)[0:2]
                xRow2Left += xSampled[1:]
                d1Row2Left += [dSampled[1]] + [findDerivativeBetweenPoints(xSampled[1], xSampled[2])]

                for n in range(1, len(xLoopsLeft)):
                    xRow2Left.append(xLoopsLeft[n][loopIdx])
                    d1Row2Left.append(findDerivativeBetweenPoints(xLoopsLeft[n-1][loopIdx], xLoopsLeft[n][loopIdx]))

                xRow2Left.append(xUp[0][-1])
                d1Row2Left.append(findDerivativeBetweenPoints(xLoopsLeft[-1][1], xUp[0][-1]))

        # Smooth derivatives from triple point to 6 point junction
        # Start from GC at upstream bifurcation ring to annulus to 6 point junction ring on right then left
        xLoopTripleTo6Pt = []
        dLoopTripleTo6Pt = []

        xLoopTripleTo6Pt += xBifurcationRings[0][0:int(len(xBifurcationRings[0]) * 0.5) + 1]
        for n2 in range(elementsAroundQuarterEso - 1):
            xLoopTripleTo6Pt.append(xAlongAround[n2][int(len(xAlongAround[n2]) * 0.5)])
            junctionIdx = n2 + 1
        xLoopTripleTo6Pt += xAlongAround[junctionIdx][int(len(xAlongAround[junctionIdx]) * 0.5):] + \
                            xAlongAround[junctionIdx][0: int(len(xAlongAround[junctionIdx]) * 0.5 + 1)]
        for n2 in range(elementsAroundQuarterEso - 1): # Note order here - going upstream
            idx = junctionIdx - 1 - n2
            xLoopTripleTo6Pt.append(xAlongAround[idx][int(len(xAlongAround[idx]) * 0.5) + 1])
        xLoopTripleTo6Pt += xBifurcationRings[0][int(len(xBifurcationRings[0]) * 0.5 + 1):]

        for n in range(len(xLoopTripleTo6Pt)):
            d = findDerivativeBetweenPoints(xLoopTripleTo6Pt[n], xLoopTripleTo6Pt[(n+1) % len(xLoopTripleTo6Pt)])
            dLoopTripleTo6Pt.append(d)
        dSmoothLoopTripleTo6Pt = interp.smoothCubicHermiteDerivativesLoop(xLoopTripleTo6Pt, dLoopTripleTo6Pt)

        # Smooth derivatives around top loop
        # Starts from GC at downstream bifurcation ring to annulus and back
        xLoopGCTriplePt = []
        dLoopGCTriplePt = []

        xLoopGCTriplePt += xBifurcationRings[1][:int(len(xBifurcationRings[1]) * 0.5) + 1]

        for n2 in range(elementsAroundQuarterEso - 2):
            idx = -(3 + n2)
            xLoopGCTriplePt.append(xUp[idx][int(len(xUp[idx]) * 0.5)])

        xLoopGCTriplePt += [xRow2Right[-1]] + [xEsoToDuodGC[1]] + [xRow2Left[0]]

        for n2 in range(elementsAroundQuarterEso - 2):
            xLoopGCTriplePt.append(xUp[n2][int(len(xUp[n2]) * 0.5) + 1])

        xLoopGCTriplePt += xBifurcationRings[1][int(len(xBifurcationRings[1]) * 0.5) + 1:]

        for n in range(len(xLoopGCTriplePt)):
            d = findDerivativeBetweenPoints(xLoopGCTriplePt[n], xLoopGCTriplePt[(n+1) % len(xLoopGCTriplePt)])
            dLoopGCTriplePt.append(d)
        dSmoothLoopGCTriplePt = interp.smoothCubicHermiteDerivativesLoop(xLoopGCTriplePt, dLoopGCTriplePt)

        # Assemble nodes and d1
        xOuter = []
        d1Outer = []
        countUp = 0
        countDown = 0

        for n2 in range(elementsCountAlong + 1):
            xAround = []
            d1Around = []
            if n2 == 0:
                for i in range(elementsAroundHalfDuod - 2):
                    xAround.append(xEsoToDuodGC[i + 1])
                    d1Around.append(d2EsoToDuodGC[i + 1])

            elif n2 == 1:
                xAround = [xEsoToDuodGC[i + n2 + 1]] + xRow2Right[1:] + xRow2Left[:-1]
                d1Around = [d2EsoToDuodGC[i + n2 + 1]] + d1Row2Right[1:] + d1Row2Left[:-1]

            elif n2 > 1 and n2 < elementsAroundQuarterEso + 2:
                xAround = xUp[countUp]
                if n2 < elementsAroundQuarterEso: # upstream of triple pt
                    d1Around = d1Up[countUp]
                    d1Around = smoothD1Around(xAround, d1Around)

                elif n2 == elementsAroundQuarterEso: # upstream bifurcation
                    # take smoothed d1 from dSmoothTripleTo6Pt
                    d1Around = dSmoothLoopTripleTo6Pt[: int(len(xBifurcationRings[0]) * 0.5) + 1] + \
                               dSmoothLoopTripleTo6Pt[-int(len(xBifurcationRings[0]) * 0.5) : ]

                elif n2 > elementsAroundQuarterEso: # downstream bifurcation
                    # take smoothed d1 from dSmoothGCToTriplePt
                    d1Around = dSmoothLoopGCTriplePt[: int(len(xBifurcationRings[1]) * 0.5) + 1] + \
                               dSmoothLoopGCTriplePt[-int(len(xBifurcationRings[1]) * 0.5) : ]
                countUp += 1

            elif n2 > elementsAroundQuarterEso + 1:
                xAround = xAlongAround[countDown]
                d1Around = d1AlongAround[countDown]

                if n2 < elementsAroundHalfEso + 1:
                    d1Around = smoothD1Around(xAround, d1Around)

                elif n2 == elementsAroundHalfEso + 1: # 6 point junction ring
                    # take smoothed d1 from dSmoothedTripleTo6Pt
                    startRightIdx = int(len(xBifurcationRings[0]) * 0.5 + elementsAroundQuarterEso + len(xAlongAround[junctionIdx]) * 0.5)
                    endRightIdx = startRightIdx + int(len(xAlongAround[junctionIdx]) * 0.5) + 1
                    startLeftIdx = startRightIdx - int(len(xAlongAround[junctionIdx]) * 0.5) + 1
                    d1Around = dSmoothLoopTripleTo6Pt[startRightIdx: endRightIdx] + \
                               dSmoothLoopTripleTo6Pt[startLeftIdx : startRightIdx]
                countDown += 1

            xOuter.append(xAround)
            d1Outer.append(d1Around)

        # Calculate d2
        xRegularLoops = []
        d2RegularLoops = []
        d2RegularOrderedLoops = []

        for n1 in range(elementsAroundHalfDuod - 2):
            xRegularLoop = []
            d1RegularRightLoop = []
            d2RegularLoop = []
            for n2 in range(elementsCountAlong):
                idx = -(1 + n2)
                xRegularLoop.append(xOuter[idx][int(len(xOuter[idx]) * 0.5 - 1 - n1)])
                d1RegularRightLoop.append(d1Outer[idx][int(len(xOuter[idx]) * 0.5 - 1 - n1)])
            xRegularLoop.append(xEsoToDuodGC[n1 + 2])
            for n2 in range(elementsCountAlong):
                xRegularLoop.append(xOuter[n2 + 1][int(len(xOuter[n2 + 1]) * 0.5 + n1 + (1 if n2 >= elementsAroundHalfEso else 2))])

            for n in range(len(xRegularLoop) - 1):
                d = findDerivativeBetweenPoints(xRegularLoop[n], xRegularLoop[n+1])
                d2RegularLoop.append(d)
            d2RegularLoop.append(d)

            d2SmoothRegularLoop = interp.smoothCubicHermiteDerivativesLine(xRegularLoop, d2RegularLoop)
            d2SmoothRegularOrderedLoop = copy.deepcopy(d2SmoothRegularLoop)

            # Switch direction on right side
            for n2 in range(elementsCountAlong):
                rotAxis = vector.normalise(vector.crossproduct3(vector.normalise(d1RegularRightLoop[n2]), d2SmoothRegularLoop[n2]))
                rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
                d = d2SmoothRegularLoop[n2]
                d2SmoothRegularLoop[n2] = [rotFrame[j][0] * d[0] + rotFrame[j][1] * d[1] + rotFrame[j][2] * d[2] for j in range(3)]
            xRegularLoops.append(xRegularLoop)
            d2RegularLoops.append(d2SmoothRegularLoop)
            d2RegularOrderedLoops.append(d2SmoothRegularOrderedLoop)

        # Smooth d2 along row 2
        xLoop2Right = []
        d1Loop2Right = []
        d2Loop2Right = []

        for n2 in range(len(xAlongAround) + len(xUp) - 1):
            idx = -(1 + n2)
            xLoop2Right.append(xOuter[idx][1])
            d1Loop2Right.append(d1Outer[idx][1])
        xLoop2Right += xRow2Right
        d1Loop2Right += d1Row2Right

        for n in range(len(xLoop2Right) - 1):
            d = findDerivativeBetweenPoints(xLoop2Right[n], xLoop2Right[n + 1])
            d2Loop2Right.append(d)
        d2Loop2Right.append(d1Row2Right[-1])
        d2Loop2Right = interp.smoothCubicHermiteDerivativesLine(xLoop2Right, d2Loop2Right, fixEndDirection=True)

        # Switch direction of d2 for downstream nodes
        for n2 in range(len(xAlongAround) + len(xUp)):
            rotAxis = vector.normalise(
                vector.crossproduct3(vector.normalise(d1Loop2Right[n2]), d2Loop2Right[n2]))
            rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
            d = d2Loop2Right[n2]
            d2Loop2Right[n2] = [rotFrame[j][0] * d[0] + rotFrame[j][1] * d[1] + rotFrame[j][2] * d[2] for j in range(3)]
            idxSwitchToD1 = n2

        # Left
        xLoop2Left = []
        d2Loop2Left = []
        xLoop2Left += xRow2Left
        for n2 in range(3, len(xOuter)):
            xLoop2Left.append(xOuter[n2][-1])

        d2Loop2Left.append(d1Row2Left[0])
        for n in range(1, len(xLoop2Left) - 1):
            d = findDerivativeBetweenPoints(xLoop2Left[n], xLoop2Left[n + 1])
            d2Loop2Left.append(d)
        d2Loop2Left.append(d)

        d2Loop2Left = interp.smoothCubicHermiteDerivativesLine(xLoop2Left, d2Loop2Left, fixStartDirection=True)

        # Smooth lower curvature
        xLC = []
        d2LC = []
        for n2 in range(elementsAroundHalfEso + 1, elementsCountAlong + 1):
            xLC.append(xOuter[n2][int(len(xOuter[n2]) * 0.5)])

        for n in range(len(xLC) - 1):
            d = findDerivativeBetweenPoints(xLC[n], xLC[n + 1])
            d2LC.append(d)
        d2LC.append(d)

        d2LC = interp.smoothCubicHermiteDerivativesLine(xLC, d2LC, fixStartDirection=True)

        # Smooth greater curvature
        d2GC = []
        for n in range(len(xEsoToDuodGC) - 1):
            d = findDerivativeBetweenPoints(xEsoToDuodGC[n], xEsoToDuodGC[n + 1])
            d2GC.append(d)
        d2GC.append(d)
        d2GC = interp.smoothCubicHermiteDerivativesLine(xEsoToDuodGC, d2GC, fixStartDirection=True)

        # Update d1 for upstream nodes
        for n1 in range(1, len(xRow2Right)):
            d1Outer[1][n1] = d2Loop2Right[idxSwitchToD1 + n1]
        for n1 in range(1, len(xRow2Left)):
            d1Outer[1][int(len(d1Outer[1]) * 0.5) + n1] = d2Loop2Left[n1 - 1]

        # Assemble d2
        d2Outer = []
        for n2 in range(elementsCountAlong + 1):
            d2Around = []
            if n2 == 0:
                d2Around.append(dSmoothLoopGCTriplePt[int(len(dSmoothLoopGCTriplePt) * 0.5)])
                for n1 in range(len(xOuter[0]) - 1):
                    d2Around.append(d2RegularLoops[n1][int(len(xRegularLoops[n1]) * 0.5)])
                    nextIdx = n1 + 1

            elif n2 == 1: # Row 2
                d2Around.append(d2RegularLoops[nextIdx][int(len(xRegularLoops[nextIdx]) * 0.5)])

                for n1 in range(nextIdx, -1, -1):
                    d2Around.append(d2RegularLoops[n1][int(len(d2RegularLoops[n1]) * 0.5) - n2])

                # right point on annulus
                d2 = dSmoothLoopGCTriplePt[int(len(xLoopGCTriplePt) * 0.5) - n2]
                rotAxis = vector.normalise(vector.crossproduct3(vector.normalise(d1Outer[n2][int(len(d1Outer[n2]) * 0.5)]), vector.normalise(d2)))
                rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
                d2Around.append([rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)])

                # left point on annulus
                d2Around.append(dSmoothLoopGCTriplePt[int(len(xLoopGCTriplePt) * 0.5) + n2])

                for n1 in range(nextIdx + 1):
                    d2Around.append(d2RegularLoops[n1][int(len(d2RegularLoops[n1]) * 0.5) + n2])

            elif n2 > 1 and n2 < elementsAroundQuarterEso + 2:
                # GC before triple point & triple point
                d2Around.append(d2GC[len(xOuter[0]) + n2])

                # Row 2 right
                d2Around.append(d2Loop2Right[-(len(xOuter[0]) + n2)])

                # Regular up right
                for n1 in range(nextIdx, -1, -1):
                    d2Around.append(d2RegularLoops[n1][int(len(d2RegularLoops[n1]) * 0.5) - n2])

                # Annulus right
                d2 = dSmoothLoopGCTriplePt[int(len(xLoopGCTriplePt) * 0.5) - n2 + (1 if n2 > elementsAroundQuarterEso else 0)]
                if n2 <= elementsAroundQuarterEso: # Rotate to point towards duodenum
                    rotAxis = vector.normalise(vector.crossproduct3(vector.normalise(d1Outer[n2][int(len(d1Outer[n2]) * 0.5)]),
                                             vector.normalise(d2)))
                    rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
                    d2Around.append(
                        [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)])
                else:
                    d2Around.append(d2) # just take d2 as-is cos we are going to remove this point later

                # Annulus left
                d2Around.append(dSmoothLoopGCTriplePt[int(len(xLoopGCTriplePt) * 0.5) + n2 - (1 if n2 > elementsAroundQuarterEso else 0)])

                # Regular down left
                for n1 in range(nextIdx + 1):
                    d2Around.append(d2RegularLoops[n1][int(len(d2RegularLoops[n1]) * 0.5) + n2])

                # Row 2 left
                d2Around.append(d2Loop2Left[len(xOuter[0]) + n2 - 1])

            elif n2 > elementsAroundQuarterEso + 1:
                # GC downstream of triple point
                d2Around.append(d2GC[len(xOuter[0]) + n2])

                # Row 2 right
                d2Around.append(d2Loop2Right[-(len(xOuter[0]) + n2)])

                # Regular up right
                for n1 in range(nextIdx, -1, -1):
                    d2Around.append(d2RegularLoops[n1][int(len(d2RegularLoops[n1]) * 0.5) - n2])

                if n2 <= elementsAroundHalfEso + 1:
                    # Annulus right between triple and 6 pt
                    idx = int(len(xBifurcationRings[0]) * 0.5 + n2 - elementsAroundQuarterEso - 1)
                    if n2 == elementsAroundHalfEso + 1:
                        d1 = dSmoothLoopTripleTo6Pt[idx]
                        d1Outer[n2][int(len(d1Outer[n2]) * 0.5)] = d1
                    else:
                        d2Around.append(dSmoothLoopTripleTo6Pt[idx])

                    # Annulus left - Rotated to point towards duodenum
                    d2 = dSmoothLoopTripleTo6Pt[-idx]
                    if n2 < elementsAroundHalfEso + 1:
                        rotAxis = vector.normalise(
                            vector.crossproduct3(vector.normalise(d1Outer[n2][int(len(d1Outer[n2]) * 0.5 + 1)]),
                                                 vector.normalise(d2)))
                    else: # use d2 on previous overlapping point to rotate
                        rotAxis = vector.normalise(
                            vector.crossproduct3(vector.normalise(d1), vector.normalise(d2)))
                    rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
                    d2Around.append(
                        [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in
                         range(3)])

                elif n2 > elementsAroundHalfEso + 1:
                    # LC - beyond 6 pt junction
                    d2Around.append(d2LC[n2 - (elementsAroundHalfEso + 1)])

                # Regular down left
                for n1 in range(nextIdx + 1):
                    d2Around.append(d2RegularLoops[n1][int(len(d2RegularLoops[n1]) * 0.5) + n2])

                # Row 2 left
                d2Around.append(d2Loop2Left[len(xOuter[0]) + n2 - 1])
            d2Outer.append(d2Around)

        # remove triple point on both sides from downstream ring
        n2Idx = elementsAroundQuarterEso + 1
        n1Idx = int(len(xOuter[n2Idx]) * 0.5)
        del xOuter[n2Idx][n1Idx: n1Idx + 2], d1Outer[n2Idx][n1Idx: n1Idx + 2], d2Outer[n2Idx][n1Idx: n1Idx + 2]

        d3UnitOuter = []
        for n2 in range(elementsCountAlong + 1):
            d3Around = []
            for n1 in range(len(xOuter[n2])):
                d3Around.append(vector.normalise(
                    vector.crossproduct3(vector.normalise(d1Outer[n2][n1]), vector.normalise(d2Outer[n2][n1]))))
            d3UnitOuter.append(d3Around)

        # Calculate curvatures
        # Curvature along GC
        xGC = []
        dGC = []
        norms = []
        for n1 in range(len(xOuter[0])):
            xGC.append(xOuter[0][n1])
            dGC.append(d1Outer[0][n1])
            norms.append(d3UnitOuter[0][n1])
        for n2 in range(1, elementsCountAlong + 1):
            xGC.append(xOuter[n2][0])
            dGC.append(d1Outer[n2][0] if n2 == 1 else d2Outer[n2][0])
            norms.append(d3UnitOuter[n2][0])
        curvatureAlongGC = findCurvatureAlongLine(xGC, dGC, norms)  # 1st len(xOuter[0]) + 1 are for d1, the rest for d2

        # Curvature along rows adjacent to GC - calculate with left and use for right as well
        norms = []
        for n in range(int(len(xOuter[1]) * 0.5)):  # d1
            norms.append(d3UnitOuter[1][n + int(len(xOuter[1]) * 0.5) + 1])
        for n2 in range(2, elementsCountAlong + 1):  # d2
            norms.append(d3UnitOuter[n2][-1])
        curvatureAlong2Left = findCurvatureAlongLine(xLoop2Left, d2Loop2Left, norms)

        # Curvature along LC
        norms = []
        for n in range(elementsAroundHalfEso + 1, elementsCountAlong + 1):
            norms.append(d3UnitOuter[n][int(len(xOuter[n]) * 0.5)])
        curvatureAlongLC = findCurvatureAlongLine(xLC[1:], d2LC[1:], norms)

        # Curvature along path from triple point to 6 point junction
        norms = []
        idxToAnnulus = elementsAroundHalfDuod + 1
        norms += d3UnitOuter[elementsAroundQuarterEso][:idxToAnnulus]

        for n2 in range(elementsAroundQuarterEso):
            idx = elementsAroundQuarterEso + 2 + n2
            if idx < elementsAroundHalfEso + 1:
                norms.append(d3UnitOuter[idx][idxToAnnulus - 1])
        norms += d3UnitOuter[elementsAroundHalfEso + 1][idxToAnnulus - 1:] + \
                 d3UnitOuter[elementsAroundHalfEso + 1][: idxToAnnulus]

        for n2 in range(elementsAroundQuarterEso - 1):
            idx = elementsAroundHalfEso - n2
            norms.append(d3UnitOuter[idx][idxToAnnulus])
        norms += d3UnitOuter[elementsAroundQuarterEso][idxToAnnulus:]
        curvatureLoopTripleTo6Pt = findCurvatureAroundLoop(xLoopTripleTo6Pt, dSmoothLoopTripleTo6Pt, norms)

        # Curvature along path from GC to triple point
        norms = []
        norms += d3UnitOuter[elementsAroundQuarterEso + 1][:idxToAnnulus - 1]
        for n2 in range(elementsAroundQuarterEso):
            idx = elementsAroundQuarterEso - n2
            norms.append(d3UnitOuter[idx][int(len(xOuter[idx]) * 0.5)])
        norms.append(d3UnitOuter[0][0])
        for n2 in range(1, elementsAroundQuarterEso + 1):
            norms.append(d3UnitOuter[n2][int(len(xOuter[n2]) * 0.5) + 1])
        norms += d3UnitOuter[elementsAroundQuarterEso + 1][idxToAnnulus - 1:]
        curvatureLoopGCTriplePt = findCurvatureAroundLoop(xLoopGCTriplePt, dSmoothLoopGCTriplePt, norms)

        # Curvature around regular loops
        curvatureRegularLoops = []
        for n1 in range(elementsAroundHalfDuod - 2):
            norms = []
            for n2 in range(elementsCountAlong):
                idx = -(1 + n2)
                norms.append(d3UnitOuter[idx][int(len(xOuter[idx]) * 0.5 - 1 - n1)])
            if n1 < elementsAroundHalfDuod - 3:
                norms.append(d3UnitOuter[0][n1 + 1])
            else:
                norms.append(d3UnitOuter[idx][0])
            for n2 in range(elementsCountAlong):
                norms.append(d3UnitOuter[n2 + 1][int(
                    len(xOuter[n2 + 1]) * 0.5 + n1 + (1 if n2 >= elementsAroundHalfEso else 2))])
            curvatureLoop = findCurvatureAlongLine(xRegularLoops[n1], d2RegularOrderedLoops[n1], norms)
            curvatureRegularLoops.append(curvatureLoop)

        # Assemble curvatures
        d1Curvature = []
        d2Curvature = []
        for n2 in range(elementsCountAlong + 1):
            d1CurvatureAround = []
            d2CurvatureAround = []
            if n2 == 0: # GC
                for i in range(elementsAroundHalfDuod - 2):
                    d1CurvatureAround.append(curvatureAlongGC[i])
                d2CurvatureAround.append(curvatureLoopGCTriplePt[int(len(curvatureLoopGCTriplePt) * 0.5)])
                for n1 in range(len(xOuter[0]) - 1):
                    d2CurvatureAround.append(curvatureRegularLoops[n1][int(len(curvatureRegularLoops[n1]) * 0.5)])
                    nextIdx = n1 + 1
            elif n2 == 1: #Row 2
                d1CurvatureAround.append(curvatureAlongGC[i + n2])
                for n in range(int(len(xOuter[1]) * 0.5) - 1, -1, -1):
                    d1CurvatureAround.append(curvatureAlong2Left[n])
                d1CurvatureAround += curvatureAlong2Left[:int(len(xOuter[1]) * 0.5)]
                d2CurvatureAround.append(curvatureRegularLoops[nextIdx][int(len(curvatureRegularLoops[nextIdx]) * 0.5)])

                for n1 in range(nextIdx, -1, -1):
                    d2CurvatureAround.append(curvatureRegularLoops[n1][int(len(curvatureRegularLoops[n1]) * 0.5) - n2])
                # right point on annulus
                d2CurvatureAround.append(curvatureLoopGCTriplePt[int(len(curvatureLoopGCTriplePt) * 0.5) - n2])
                # left point on annulus
                d2CurvatureAround.append(curvatureLoopGCTriplePt[int(len(curvatureLoopGCTriplePt) * 0.5) + n2])
                for n1 in range(nextIdx + 1):
                    d2CurvatureAround.append(curvatureRegularLoops[n1][int(len(curvatureRegularLoops[n1]) * 0.5) + n2])

            elif n2 > 1 and n2 < elementsAroundQuarterEso + 2: # Before triple pt & triple point
                xAround = xOuter[n2]
                if n2 < elementsAroundQuarterEso:  # upstream of triple pt
                    d1Around = d1Outer[n2]
                    normsAround = d3UnitOuter[n2]
                    d1CurvatureAround = findD1CurvatureAround(xAround, d1Around, normsAround)

                elif n2 == elementsAroundQuarterEso:  # upstream bifurcation
                    # take smoothed d1 from dSmoothTripleTo6Pt
                    d1CurvatureAround = curvatureLoopTripleTo6Pt[: int(len(xBifurcationRings[0]) * 0.5) + 1] + \
                                        curvatureLoopTripleTo6Pt[-int(len(xBifurcationRings[0]) * 0.5):]

                elif n2 > elementsAroundQuarterEso:  # downstream bifurcation
                    # take smoothed d1 from dSmoothGCToTriplePt
                    d1CurvatureAround = curvatureLoopGCTriplePt[: int(len(xBifurcationRings[1]) * 0.5) + 1] + \
                                        curvatureLoopGCTriplePt[-int(len(xBifurcationRings[1]) * 0.5):]

                # GC
                d2CurvatureAround.append(curvatureAlongGC[len(xOuter[0]) + n2 - 1])
                # Row 2 right
                d2CurvatureAround.append(curvatureAlong2Left[len(xOuter[0]) + n2 - 1])
                # Regular up right
                for n1 in range(nextIdx, -1, -1):
                    d2CurvatureAround.append(curvatureRegularLoops[n1][int(len(curvatureRegularLoops[n1]) * 0.5) - n2])
                # Annulus right
                d2CurvatureAround.append(curvatureLoopGCTriplePt[
                    int(len(curvatureLoopGCTriplePt) * 0.5) - n2 + (1 if n2 > elementsAroundQuarterEso else 0)])
                # Annulus left
                d2CurvatureAround.append(curvatureLoopGCTriplePt[int(len(curvatureLoopGCTriplePt) * 0.5) + n2 - (
                    1 if n2 > elementsAroundQuarterEso else 0)])
                # Regular down left
                for n1 in range(nextIdx + 1):
                    d2CurvatureAround.append(curvatureRegularLoops[n1][int(len(curvatureRegularLoops[n1]) * 0.5) + n2])
                # Row 2 left
                d2CurvatureAround.append(curvatureAlong2Left[len(xOuter[0]) + n2 - 1])

            elif n2 > elementsAroundQuarterEso + 1: # Downstream of triple point
                xAround = xOuter[n2]
                d1Around = d1Outer[n2]
                normsAround = d3UnitOuter[n2]

                if n2 < elementsAroundHalfEso + 1:
                    d1CurvatureAround = findD1CurvatureAround(xAround, d1Around, normsAround)

                elif n2 == elementsAroundHalfEso + 1:  # 6 point junction ring
                    # take smoothed d1 from dSmoothedTripleTo6Pt
                    startRightIdx = int(len(xBifurcationRings[0]) * 0.5 + elementsAroundQuarterEso + len(
                        xAlongAround[junctionIdx]) * 0.5)
                    endRightIdx = startRightIdx + int(len(xAlongAround[junctionIdx]) * 0.5) + 1
                    startLeftIdx = startRightIdx - int(len(xAlongAround[junctionIdx]) * 0.5) + 1
                    d1CurvatureAround = curvatureLoopTripleTo6Pt[startRightIdx: endRightIdx] + \
                                        curvatureLoopTripleTo6Pt[startLeftIdx: startRightIdx]

                if n2 > elementsAroundHalfEso + 1: # closed rings beyond 6 point junction
                    xLoop = xAround[int(len(xAround) * 0.5 + 1):] + xAround[: int(len(xAround) * 0.5 + 1)]
                    d1Loop = d1Around[int(len(d1Around) * 0.5 + 1):] + d1Around[: int(len(d1Around) * 0.5 + 1)]
                    normsLoop = normsAround[int(len(normsAround) * 0.5 + 1):] + normsAround[ : int(len(normsAround) * 0.5 + 1)]
                    curvature = findCurvatureAroundLoop(xLoop, d1Loop, normsLoop)
                    # Rearrange to correct order
                    d1CurvatureAround = curvature[int(len(xLoop) * 0.5) - 1:] + curvature[: int(len(xAround) * 0.5) - 1]

                # GC
                d2CurvatureAround.append(curvatureAlongGC[len(xOuter[0]) + n2 - 1])
                # Row 2 right
                d2CurvatureAround.append(curvatureAlong2Left[len(xOuter[0]) + n2 - 1])
                # Regular up right
                for n1 in range(nextIdx, -1, -1):
                    d2CurvatureAround.append(curvatureRegularLoops[n1][int(len(curvatureRegularLoops[n1]) * 0.5) - n2])
                if n2 <= elementsAroundHalfEso + 1:
                    # Annulus right between triple and 6 pt
                    idx = int(len(xBifurcationRings[0]) * 0.5 + n2 - elementsAroundQuarterEso - 1)
                    if n2 == elementsAroundHalfEso + 1:
                        d1CurvatureAround[int(len(d1Outer[n2]) * 0.5)] = curvatureLoopTripleTo6Pt[idx]
                    else:
                        d2CurvatureAround.append(curvatureLoopTripleTo6Pt[idx])
                    # Annulus left
                    d2CurvatureAround.append(curvatureLoopTripleTo6Pt[-idx])
                elif n2 > elementsAroundHalfEso + 1: # Beyond 6 pt junction
                    # LC
                    d2CurvatureAround.append(curvatureAlongLC[n2 - (elementsAroundHalfEso + 1) - 1])
                # Regular down left
                for n1 in range(nextIdx + 1):
                    d2CurvatureAround.append(curvatureRegularLoops[n1][int(len(curvatureRegularLoops[n1]) * 0.5) + n2])
                # Row 2 left
                d2CurvatureAround.append(curvatureAlong2Left[len(xOuter[0]) + n2 - 1])
            d1Curvature.append(d1CurvatureAround)
            d2Curvature.append(d2CurvatureAround)

        # Create inner nodes
        xList = []
        d1List = []
        d2List = []
        d3List = []
        nodeIdx = stomachStartNode
        idxMat = []

        for n2 in range(elementsCountAlong + 1):
            idxThroughWall = []
            for n3 in range(elementsCountThroughWall + 1):
                xi3 = 1 / elementsCountThroughWall * n3
                idxAround = []
                for n1 in range(len(xOuter[n2])):
                    # Coordinates
                    norm = d3UnitOuter[n2][n1]
                    xOut = xOuter[n2][n1]
                    xIn = [xOut[i] - norm[i] * wallThickness for i in range(3)]
                    dWall = [wallThickness * c for c in norm]
                    x = interp.interpolateCubicHermite(xIn, dWall, xOut, dWall, xi3)
                    xList.append(x)

                    # d1
                    factor = 1.0 + wallThickness * (1.0 - xi3) * d1Curvature[n2][n1]
                    d1 = [factor * c for c in d1Outer[n2][n1]]
                    d1List.append(d1)

                    # d2
                    factor = 1.0 + wallThickness * (1.0 - xi3) * d2Curvature[n2][n1]
                    d2 = [factor * c for c in d2Outer[n2][n1]]
                    d2List.append(d2)

                    # d3
                    d3 = [c * wallThickness / elementsCountThroughWall for c in norm]
                    d3List.append(d3)

                    idxAround.append(nodeIdx)
                    nodeIdx += 1
                idxThroughWall.append(idxAround)
            idxMat.append(idxThroughWall)

        nodeIdxGC = []
        nodesFlipD2 = []
        for n2 in range(len(idxMat)):
            for n3 in range(len(idxMat[n2])):
                if n2 == 0:
                    nodeIdxGC += idxMat[n2][n3]
                    nodesFlipD2 += idxMat[n2][n3]
                else:
                    nodeIdxGC.append(idxMat[n2][n3][0])

        nodeIdxLC = []
        for n2 in range(elementsCountAlong - elementsAlongCardiaToDuod, elementsCountAlong + 1):
            for n3 in range(len(idxMat[n2])):
                nodeIdxLC.append(idxMat[n2][n3][elementsAroundHalfDuod])

        for n2 in range(len(xList)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xList[n2])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1List[n2])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2List[n2])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3List[n2])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier += 1

        # Create element
        elementIdxMat = []
        n = 0
        for n2 in range(elementsAlongEsophagus):
            elementIdxThroughWall = []
            for n3 in range(elementsThroughEsophagusWall):
                elementIdxAround = []
                for n1 in range(elementsCountAroundEso):
                    n += 1
                    elementIdxAround.append(n)
                elementIdxThroughWall.append(elementIdxAround)
            elementIdxMat.append(elementIdxThroughWall)

        if useCubicHermiteThroughWall:
            eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        else:
            eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eftStandard = eftfactory.createEftBasic()

        elementtemplateStandard = mesh.createElementtemplate()
        elementtemplateStandard.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateStandard.defineField(coordinates, -1, eftStandard)

        elementtemplateX = mesh.createElementtemplate()
        elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        for e2 in range(elementsCountAlong):
            startNode = stomachStartNode
            for e in range(e2):
                startNode += len(xOuter[e]) * (elementsCountThroughWall + 1)
            elementsCountAround1 = len(xOuter[e2])
            elementsAroundThroughWall = elementsCountAround1 * (elementsCountThroughWall + 1)
            elementsCountAround2 = len(xOuter[e2 + 1])

            # Row 1
            if e2 == 0:
                elementIdxThroughWall = []
                for e3 in range(elementsCountThroughWall):
                    elementIdxAround = []
                    for e1 in range(int(elementsCountAround1) * 2 + 1):
                        if e1 != elementsCountAround1:
                            scaleFactors = []
                            eft1 = eftStandard
                            elementtemplate1 = elementtemplateStandard
                            if e1 < elementsCountAround1:
                                if e1 == 0:
                                    bni11 = startNode + elementsAroundThroughWall + e3 * elementsCountAround2 + e1
                                    bni12 = startNode + elementsCountAround1 - e1 + e3 * elementsCountAround1 - 1
                                else:
                                    bni11 = startNode + elementsCountAround1 - e1 + e3 * elementsCountAround1
                                    bni12 = bni11 - 1
                                bni21 = startNode + elementsAroundThroughWall + 1 + e1 + e3 * elementsCountAround2
                                bni22 = bni21 + 1
                                nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                                   bni11 + (elementsCountAround2 if e1 == 0 else elementsCountAround1),
                                                   bni12 + elementsCountAround1,
                                                   bni21 + elementsCountAround2, bni22 + elementsCountAround2]
                                eft1 = eftfactory.createEftNoCrossDerivatives()
                                scaleFactors = [-1.0]
                                setEftScaleFactorIds(eft1, [1], [])
                                scaleEftNodeValueLabels(eft1, [1, 2, 5, 6],
                                                        [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D2_DS1DS2,
                                                         Node.VALUE_LABEL_D2_DS1DS3,
                                                         Node.VALUE_LABEL_D3_DS1DS2DS3], [1])
                                scaleEftNodeValueLabels(eft1, [1, 2, 5, 6],
                                                        [Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                                                         Node.VALUE_LABEL_D2_DS2DS3,
                                                         Node.VALUE_LABEL_D3_DS1DS2DS3], [1])
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX

                            elif e1 > elementsCountAround1:
                                if e1 < elementsCountAround1 * 2:
                                    bni11 = startNode + e1 - elementsCountAround1 - 1 + elementsCountAround1 * e3
                                    bni12 = bni11 + 1
                                else:
                                    bni11 = startNode + elementsCountAround1 + e3 * elementsCountAround1 - 1
                                    bni12 = startNode + elementsAroundThroughWall + e3 * elementsCountAround2
                                bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3 + 1
                                bni22 = bni21 + 1
                                nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                                   bni11 + elementsCountAround1,
                                                   bni12 + (elementsCountAround1 if e1 < elementsCountAround1 * 2 else elementsCountAround2),
                                                   bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                            element = mesh.createElement(elementIdentifier, elementtemplate1)
                            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                            if scaleFactors:
                                result3 = element.setScaleFactors(eft1, scaleFactors)
                            elementIdxAround.append(elementIdentifier)
                            elementIdentifier += 1
                            annotationGroups = annotationGroupsAlong[e2]
                            if annotationGroups:
                                allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                                for annotationGroup in annotationGroups:
                                    meshGroup = annotationGroup.getMeshGroup(mesh)
                                    meshGroup.addElement(element)
                    elementIdxThroughWall.append(elementIdxAround)
                elementIdxMat.append(elementIdxThroughWall)

            # Row 2
            elif e2 == 1:
                elementIdxThroughWall = []
                for e3 in range(elementsCountThroughWall):
                    elementIdxAround = []
                    for e1 in range(elementsCountAround1 + 2):
                        if e1 != int(elementsCountAround1 * 0.5 + 1):
                            scaleFactors = []
                            eft1 = eftStandard
                            elementtemplate1 = elementtemplateStandard
                            if e1 < 2:
                                bni11 = startNode + e3 * elementsCountAround1 + e1
                                bni12 = startNode + e3 * elementsCountAround1 + (e1 + 1)
                                bni21 = startNode + elementsAroundThroughWall + elementsCountAround2 * e3 + e1
                                bni22 = startNode + elementsAroundThroughWall + elementsCountAround2 * e3 + (e1 + 1)
                                if e1 == 0:  # Remap derivatives of element adjacent to GC
                                    scaleFactors = [-1.0]
                                    nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                                       bni11 + elementsCountAround1,
                                                       bni12 + elementsCountAround1,
                                                       bni21 + elementsCountAround2,
                                                       bni22 + elementsCountAround2]
                                    eft1 = eftfactory.createEftNoCrossDerivatives()
                                    setEftScaleFactorIds(eft1, [1], [])
                                    remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,[(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                                elif e1 == 1:  # Bottom right wedge
                                    scaleFactors = [-1.0]
                                    nodeIdentifiers = [bni11, bni21, bni22,
                                                       bni11 + elementsCountAround1,
                                                       bni21 + elementsCountAround2,
                                                       bni22 + elementsCountAround2]
                                    eft1 = eftfactory.createEftWedgeCollapseXi1Quadrant([1, 5])
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX

                            elif e1 > 1 and e1 < elementsCountAround1:
                                bni11 = startNode + e3 * elementsCountAround1 + e1 - 1
                                bni12 = startNode + e3 * elementsCountAround1 + e1 % elementsCountAround1
                                bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                                bni22 = startNode + elementsAroundThroughWall + (e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3
                                nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                                   bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                                   bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                            elif e1 >= elementsCountAround1:
                                bni11 = startNode + e3 * elementsCountAround1 + e1 - 2
                                bni12 = startNode + e3 * elementsCountAround1 + (e1 - 1) % elementsCountAround1
                                bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                                bni22 = startNode + elementsAroundThroughWall + (e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3
                                if e1 == elementsCountAround1:  # Bottom left wedge
                                    nodeIdentifiers = [bni12, bni21, bni22,
                                                       bni12 + elementsCountAround1,
                                                       bni21 + elementsCountAround2,
                                                       bni22 + elementsCountAround2]
                                    eft1 = eftfactory.createEftWedgeCollapseXi1Quadrant([2, 6])
                                elif e1 == elementsCountAround1 + 1:  # Remap derivatives of element adjacent to GC
                                    scaleFactors = [-1.0]
                                    nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                                       bni11 + elementsCountAround1,
                                                       bni12 + elementsCountAround1,
                                                       bni21 + elementsCountAround2,
                                                       bni22 + elementsCountAround2]
                                    eft1 = eftfactory.createEftNoCrossDerivatives()
                                    setEftScaleFactorIds(eft1, [1], [])
                                    remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX

                            element = mesh.createElement(elementIdentifier, elementtemplate1)
                            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                            if scaleFactors:
                                result3 = element.setScaleFactors(eft1, scaleFactors)
                            elementIdxAround.append(elementIdentifier)
                            elementIdentifier += 1
                            annotationGroups = annotationGroupsAlong[e2]
                            if annotationGroups:
                                allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                                for annotationGroup in annotationGroups:
                                    meshGroup = annotationGroup.getMeshGroup(mesh)
                                    meshGroup.addElement(element)
                    elementIdxThroughWall.append(elementIdxAround)
                elementIdxMat.append(elementIdxThroughWall)

            # Additional elements between second and upstream bifurcation ring
            elif e2 > 1 and e2 < elementsAroundQuarterEso:
                elementIdxThroughWall = []
                for e3 in range(elementsCountThroughWall):
                    elementIdxAround = []
                    for e1 in range(elementsCountAround1):
                        if e1 != int(elementsCountAround1 * 0.5):
                            scaleFactors = []
                            eft1 = eftStandard
                            elementtemplate1 = elementtemplateStandard
                            bni11 = startNode + e3 * elementsCountAround1 + e1
                            bni12 = startNode + e3 * elementsCountAround1 + (e1 + 1) % elementsCountAround1
                            bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                            bni22 = startNode + elementsAroundThroughWall + (e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3
                            nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                               bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                               bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                            element = mesh.createElement(elementIdentifier, elementtemplate1)
                            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                            if scaleFactors:
                                result3 = element.setScaleFactors(eft1, scaleFactors)
                            elementIdxAround.append(elementIdentifier)
                            elementIdentifier += 1
                            annotationGroups = annotationGroupsAlong[e2]
                            if annotationGroups:
                                allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                                for annotationGroup in annotationGroups:
                                    meshGroup = annotationGroup.getMeshGroup(mesh)
                                    meshGroup.addElement(element)
                    elementIdxThroughWall.append(elementIdxAround)
                elementIdxMat.append(elementIdxThroughWall)

            # Upstream bifurcation
            elif e2 == elementsAroundQuarterEso:
                elementIdxThroughWall = []
                for e3 in range(elementsCountThroughWall):
                    elementIdxAround = []
                    for e1 in range(elementsCountAround1):
                        if e1 != int(elementsCountAround1 * 0.5):
                            scaleFactors = []
                            eft1 = eftStandard
                            elementtemplate1 = elementtemplateStandard
                            bni11 = startNode + e3 * elementsCountAround1 + e1
                            bni12 = startNode + e3 * elementsCountAround1 + (e1 + 1) % elementsCountAround1
                            bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                            bni22 = startNode + elementsAroundThroughWall + (e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3

                            if e1 < int(elementsCountAround1 * 0.5) - 1:
                                nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                                   bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                                   bni21 + elementsCountAround2, bni22 + elementsCountAround2]
                            elif e1 == int(elementsCountAround1 * 0.5) - 1:  # right wedge
                                nodeIdentifiers = [bni11, bni12, bni21,
                                                   bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                                   bni21 + elementsCountAround2]
                                scaleFactors = [-1.0]
                                eft1 = eftfactory.createEftWedgeCollapseXi2Quadrant([4, 8])
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX

                            elif e1 == int(elementsCountAround1 * 0.5) + 1:  # left wedge
                                bni21 = bni21 - 1
                                nodeIdentifiers = [bni11, bni12, bni21,
                                                   bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                                   bni21 + elementsCountAround2]
                                eft1 = eftfactory.createEftWedgeCollapseXi2Quadrant([3, 7])
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX

                            elif e1 > int(elementsCountAround1 * 0.5) + 1:
                                bni21 = bni21 - 2
                                bni22 = startNode + elementsAroundThroughWall + (e1 - 1) % elementsCountAround2 + elementsCountAround2 * e3
                                nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                                   bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                                   bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                            element = mesh.createElement(elementIdentifier, elementtemplate1)
                            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                            if scaleFactors:
                                result3 = element.setScaleFactors(eft1, scaleFactors)
                            if e1 == 0:
                                fundusBodyJunctionElementIdentifier = elementIdentifier
                            elementIdxAround.append(elementIdentifier)
                            elementIdentifier += 1
                            annotationGroups = annotationGroupsAlong[e2]
                            if annotationGroups:
                                allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                                for annotationGroup in annotationGroups:
                                    meshGroup = annotationGroup.getMeshGroup(mesh)
                                    meshGroup.addElement(element)
                    elementIdxThroughWall.append(elementIdxAround)
                elementIdxMat.append(elementIdxThroughWall)

            # Downstream bifurcation
            elif e2 == elementsAroundQuarterEso + 1:
                elementIdxThroughWall = []
                for e3 in range(elementsCountThroughWall):
                    elementIdxAround = []
                    for e1 in range(elementsCountAround1 + 1):
                        scaleFactors = []
                        eft1 = eftStandard
                        elementtemplate1 = elementtemplateStandard
                        if e1 < int(elementsCountAround1 * 0.5) + 1:
                            bni11 = startNode + e3 * elementsCountAround1 + e1
                        elif e1 == int(elementsCountAround1 * 0.5) + 1:
                            bni11 = startNode - len(xOuter[e2 - 1]) * (elementsCountThroughWall + 1) + e3 * len(xOuter[e2 - 1]) + e1 + 1
                        elif e1 > int(elementsCountAround1 * 0.5) + 1:
                            bni11 = startNode + e3 * elementsCountAround1 + e1 - 1

                        if e1 < int(elementsCountAround1 * 0.5):
                            bni12 = startNode + e3 * elementsCountAround1 + (e1 + 1) % elementsCountAround1
                        elif e1 == int(elementsCountAround1 * 0.5):
                            bni12 = startNode - len(xOuter[e2 - 1]) * (elementsCountThroughWall + 1) + e3 * len(xOuter[e2 - 1]) + e1 + 1
                        elif e1 > int(elementsCountAround1 * 0.5):
                            bni12 = startNode + e3 * elementsCountAround1 + e1 % elementsCountAround1

                        if e1 > int(elementsCountAround1 * 0.5):
                            bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3 + 1
                            bni22 = startNode + elementsAroundThroughWall + (e1 + 2) % elementsCountAround2 + elementsCountAround2 * e3
                        else:
                            bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                            bni22 = startNode + elementsAroundThroughWall + (e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3

                        nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                           bni11 + (len(xOuter[e2 - 1]) if e1 == int(elementsCountAround1 * 0.5) + 1 else elementsCountAround1),
                                           bni12 + (len(xOuter[e2 - 1]) if e1 == int(elementsCountAround1 * 0.5) else elementsCountAround1),
                                           bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                        if e1 == int(elementsCountAround1 * 0.5):
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                        elif e1 == int(elementsCountAround1 * 0.5) + 1:
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [ ])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                        element = mesh.createElement(elementIdentifier, elementtemplate1)
                        result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                        if scaleFactors:
                            result3 = element.setScaleFactors(eft1, scaleFactors)
                        elementIdxAround.append(elementIdentifier)
                        elementIdentifier += 1
                        annotationGroups = annotationGroupsAlong[e2]
                        if annotationGroups:
                            allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                            for annotationGroup in annotationGroups:
                                meshGroup = annotationGroup.getMeshGroup(mesh)
                                meshGroup.addElement(element)
                    elementIdxThroughWall.append(elementIdxAround)
                elementIdxMat.append(elementIdxThroughWall)

            # Rows between downstream and penultimate ring
            elif elementsAroundQuarterEso + 2 <= e2 < elementsAroundHalfEso:
                elementIdxThroughWall = []
                for e3 in range(elementsCountThroughWall):
                    elementIdxAround = []
                    for e1 in range(elementsCountAround1 - 1):
                        bni11 = startNode + e3 * elementsCountAround1 + e1 + (0 if e1 < int(elementsCountAround1 * 0.5) else 1)
                        bni12 = startNode + e3 * elementsCountAround1 + (e1 + (1 if e1 < int(elementsCountAround1 * 0.5) else 2)) % elementsCountAround1
                        bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3 + (0 if e1 < int(elementsCountAround1 * 0.5) else 1)
                        bni22 = startNode + elementsAroundThroughWall + (e1 + (1 if e1 < int(elementsCountAround1 * 0.5) else 2)) % elementsCountAround2 + elementsCountAround2 * e3
                        nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                           bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                           bni21 + elementsCountAround2, bni22 + elementsCountAround2]
                        element = mesh.createElement(elementIdentifier, elementtemplateStandard)
                        result = element.setNodesByIdentifier(eftStandard, nodeIdentifiers)
                        elementIdxAround.append(elementIdentifier)
                        elementIdentifier = elementIdentifier + 1
                        annotationGroups = annotationGroupsAlong[e2]
                        if annotationGroups:
                            allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                            for annotationGroup in annotationGroups:
                                meshGroup = annotationGroup.getMeshGroup(mesh)
                                meshGroup.addElement(element)
                    elementIdxThroughWall.append(elementIdxAround)
                elementIdxMat.append(elementIdxThroughWall)

            # Penultimate row connecting to annulus and beyond
            elif elementsAroundHalfEso <= e2:
                elementIdxThroughWall = []
                for e3 in range(elementsCountThroughWall):
                    elementIdxAround = []
                    for e1 in range(elementsCountAround1 - (1 if e2 == elementsAroundHalfEso else 0)):
                        scaleFactors = []
                        eft1 = eftStandard
                        elementtemplate1 = elementtemplateStandard
                        if e2 == elementsAroundHalfEso:
                            bni11 = startNode + e3 * elementsCountAround1 + e1 + (0 if e1 < int(elementsCountAround1 * 0.5) else 1)
                            bni12 = startNode + e3 * elementsCountAround1 + (e1 + (1 if e1 < int(elementsCountAround1 * 0.5) else 2)) % elementsCountAround1
                            # Remap elements next to annulus
                            if e1 == int(elementsCountAround1 * 0.5) - 1:
                                scaleFactors = [-1.0]
                                eft1 = eftfactory.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2, ([(Node.VALUE_LABEL_D_DS1, [])]))
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX

                        else:
                            bni11 = startNode + e3 * elementsCountAround1 + e1
                            bni12 = startNode + e3 * elementsCountAround1 + (e1 + 1) % elementsCountAround1
                        bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                        bni22 = startNode + elementsAroundThroughWall + (
                                    e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3
                        nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                           bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                           bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                        if e2 == elementsAroundHalfEso + 1:
                            if e1 == int(elementsCountAround1 * 0.5) - 1:
                                scaleFactors = [-1.0]
                                eft1 = eftfactory.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2, ([(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])]))
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX

                            elif e1 == int(elementsCountAround1 * 0.5):
                                eft1 = eftfactory.createEftNoCrossDerivatives()
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2, ([(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])]))
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX

                        element = mesh.createElement(elementIdentifier, elementtemplate1)
                        result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                        if scaleFactors:
                            result3 = element.setScaleFactors(eft1, scaleFactors)
                        elementIdxAround.append(elementIdentifier)
                        elementIdentifier += 1
                        annotationGroups = annotationGroupsAlong[e2]
                        if annotationGroups:
                            allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                            for annotationGroup in annotationGroups:
                                meshGroup = annotationGroup.getMeshGroup(mesh)
                                meshGroup.addElement(element)
                    elementIdxThroughWall.append(elementIdxAround)
                elementIdxMat.append(elementIdxThroughWall)

        # Annulus
        # Assemble endPoints for annulus
        endPoints_x = [[None] * elementsCountAroundEso, [None] * elementsCountAroundEso]
        endPoints_d1 = [[None] * elementsCountAroundEso, [None] * elementsCountAroundEso]
        endPoints_d2 = [[None] * elementsCountAroundEso, [None] * elementsCountAroundEso]
        endNode_Id = [[None] * elementsCountAroundEso, [None] * elementsCountAroundEso]
        endDerivativesMap = [[None] * elementsCountAroundEso, [None] * elementsCountAroundEso]
        endProportions = []

        thicknessIdx = [0, -1]
        for nAround in range(elementsCountAroundEso):
            for n3 in range(len(thicknessIdx)):
                if nAround == 0:
                    idx = idxMat[nAround][thicknessIdx[n3]][0]
                elif nAround <= elementsAroundQuarterEso:
                    idx = idxMat[nAround][thicknessIdx[n3]][int((len(xOuter[nAround]) - 1) * 0.5)]
                elif elementsAroundQuarterEso < nAround < elementsAroundHalfEso:
                    idx = idxMat[nAround + 1][thicknessIdx[n3]][int((len(xOuter[nAround + 1]) - 1) * 0.5)]
                elif nAround == elementsAroundHalfEso:
                    idx = idxMat[nAround + 1][thicknessIdx[n3]][int(len(xOuter[nAround + 1]) * 0.5)]
                elif nAround > elementsAroundHalfEso:
                    idx = endNode_Id[n3][elementsAroundHalfEso - (nAround - elementsAroundHalfEso)] + 1

                endPoints_x[n3][nAround] = xList[idx - stomachStartNode]
                endPoints_d1[n3][nAround] = d1List[idx - stomachStartNode]
                endPoints_d2[n3][nAround] = d2List[idx - stomachStartNode]
                endNode_Id[n3][nAround] = idx

                if n3 == len(thicknessIdx) - 1: # outer layer
                    endPosition = trackSurfaceStomach.findNearestPosition(endPoints_x[n3][nAround])
                    endProportions.append(trackSurfaceStomach.getProportion(endPosition))

        for nAround in range(elementsCountAroundEso):
            if nAround == 0:
                endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((0, -1, 0), (1, 0, 0), None)
            elif nAround == elementsAroundQuarterEso:
                endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((0, 1, 0), (-1, 1, 0), None, (1, 0, 0))
            elif 0 < nAround < elementsAroundHalfEso:
                endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((0, 1, 0), (-1, 0, 0), None)
            elif nAround == elementsAroundHalfEso:
                endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((1, 0, 0), (1, 1, 0), None, (0, -1, 0))
            elif nAround == int(elementsCountAroundEso * 0.75):
                endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((1, 0, 0), (1, 1, 0), None, (0, -1, 0))
            elif elementsAroundHalfEso < nAround < elementsCountAroundEso:
                endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((0, -1, 0), (1, 0, 0), None)

        startProportions = []
        for n in range(elementsCountAroundEso):
            startProportions.append(trackSurfaceStomach.getProportion(o1_Positions[n]))

        cardiaGroup = AnnotationGroup(region, get_stomach_term("cardia of stomach"))
        cardiaMeshGroup = cardiaGroup.getMeshGroup(mesh)
        allAnnotationGroups.append(cardiaGroup)

        lastDuodenumElementIdentifier = elementIdentifier

        nextNodeIdentifier, nextElementIdentifier = createAnnulusMesh3d(
            nodes, mesh, nodeIdentifier, elementIdentifier,
            o1_x, o1_d1, o1_d2, None, o1_NodeId, None,
            endPoints_x, endPoints_d1, endPoints_d2, None, endNode_Id, endDerivativesMap,
            elementsCountRadial = elementsCountAcrossCardia, meshGroups= [stomachMeshGroup, cardiaMeshGroup],
            tracksurface=trackSurfaceStomach, startProportions = startProportions, endProportions = endProportions,
            rescaleStartDerivatives = True, rescaleEndDerivatives = True)

        elementIdxThroughWall = []
        n = lastDuodenumElementIdentifier - 1
        for n3 in range(elementsCountThroughWall):
            elementIdxAround = []
            for n1 in range(elementsCountAroundEso):
                n += 1
                elementIdxAround.append(n)
            elementIdxThroughWall.append(elementIdxAround)
        elementIdxMat.append(elementIdxThroughWall)

        nodeIdentifier = nextNodeIdentifier

        # annotation fiducial points for embedding in whole body
        GEJLCGroup = findOrCreateAnnotationGroupForTerm(allAnnotationGroups, region, get_stomach_term("gastro-esophagal junction on lesser curvature"))
        GEJLCElement = mesh.findElementByIdentifier(stomachStartElement - elementsAroundHalfEso - 1)
        GEJLCXi = [1.0, 1.0, 1.0]
        cache.setMeshLocation(GEJLCElement, GEJLCXi)
        markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
        nodeIdentifier += 1
        cache.setNode(markerPoint)
        markerName.assignString(cache, GEJLCGroup.getName())
        markerLocation.assignMeshLocation(cache, GEJLCElement, GEJLCXi)
        for group in [stomachGroup, GEJLCGroup]:
            group.getNodesetGroup(nodes).addNode(markerPoint)

        fundusBodyJunctionGroup = findOrCreateAnnotationGroupForTerm(allAnnotationGroups, region, get_stomach_term("limiting ridge on greater curvature" if limitingRidge else "junction between fundus and body on greater curvature"))
        fundusBodyJunctionElement = mesh.findElementByIdentifier(fundusBodyJunctionElementIdentifier)
        fundusBodyJunctionXi = [0.0, 0.0 if limitingRidge else 1.0, 1.0]
        cache.setMeshLocation(fundusBodyJunctionElement, fundusBodyJunctionXi)
        markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
        nodeIdentifier += 1
        cache.setNode(markerPoint)
        markerName.assignString(cache, fundusBodyJunctionGroup.getName())
        markerLocation.assignMeshLocation(cache, fundusBodyJunctionElement, fundusBodyJunctionXi)
        for group in [stomachGroup, fundusBodyJunctionGroup]:
            group.getNodesetGroup(nodes).addNode(markerPoint)

        pylorusGCGroup = findOrCreateAnnotationGroupForTerm(allAnnotationGroups, region, get_stomach_term("pylorus on greater curvature"))
        pylorusGCElement = mesh.findElementByIdentifier(lastDuodenumElementIdentifier - (elementsCountAlongGroups[-1] + 1) * elementsCountAroundDuod)
        pylorusGCXi = [0.0, 1.0, 1.0]
        cache.setMeshLocation(pylorusGCElement, pylorusGCXi)
        markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
        nodeIdentifier += 1
        cache.setNode(markerPoint)
        markerName.assignString(cache, pylorusGCGroup.getName())
        markerLocation.assignMeshLocation(cache, pylorusGCElement, pylorusGCXi)
        for group in [stomachGroup, pylorusGCGroup]:
            group.getNodesetGroup(nodes).addNode(markerPoint)

        duodenumGCGroup = findOrCreateAnnotationGroupForTerm(allAnnotationGroups, region, get_stomach_term("duodenum on greater curvature"))
        duodenumGCElement = mesh.findElementByIdentifier(lastDuodenumElementIdentifier - elementsCountAroundDuod)
        duodenumGCXi = [0.0, 1.0, 1.0]
        cache.setMeshLocation(duodenumGCElement, duodenumGCXi)
        markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
        nodeIdentifier += 1
        cache.setNode(markerPoint)
        markerName.assignString(cache, duodenumGCGroup.getName())
        markerLocation.assignMeshLocation(cache, duodenumGCElement, duodenumGCXi)
        for group in [stomachGroup, duodenumGCGroup]:
            group.getNodesetGroup(nodes).addNode(markerPoint)

        # Create annotation groups for dorsal and ventral parts of the stomach
        dorsalGroup = AnnotationGroup(region, get_stomach_term("dorsal stomach"))
        ventralGroup = AnnotationGroup(region, get_stomach_term("ventral stomach"))
        dorsalMeshGroup = dorsalGroup.getMeshGroup(mesh)
        ventralMeshGroup = ventralGroup.getMeshGroup(mesh)

        for e2 in range(len(elementIdxMat)):
            for e3 in range(len(elementIdxMat[e2])):
                for e1 in range(len(elementIdxMat[e2][e3])):
                    elementIdx = elementIdxMat[e2][e3][e1]
                    element = mesh.findElementByIdentifier(elementIdx)
                    if e1 < 0.5 * len(elementIdxMat[e2][e3]):
                        dorsalMeshGroup.addElement(element)
                    else:
                        ventralMeshGroup.addElement(element)
        allAnnotationGroups.append(dorsalGroup)
        allAnnotationGroups.append(ventralGroup)

        # Create split coordinate field
        nodesOnSplitMargin = []
        nodesOnLCMargin = []

        for n2 in range(elementsAlongEsophagus + 1):
            for n3 in range(elementsThroughEsophagusWall + 1):
                nodeIdxOnGCMargin = 1 + n2 * (elementsThroughEsophagusWall + 1 ) * elementsCountAroundEso + n3 * elementsCountAroundEso
                nodesOnSplitMargin.append(nodeIdxOnGCMargin)
                nodeIdxOnLCMargin = 1 + elementsAroundHalfEso + n2 * (elementsThroughEsophagusWall + 1 ) * elementsCountAroundEso + n3 * elementsCountAroundEso
                nodesOnSplitMargin.append(nodeIdxOnLCMargin)
                nodesOnLCMargin.append(nodeIdxOnLCMargin)
        nodesOnSplitMargin += nodeIdxGC + nodeIdxLC

        splitCoordinates = findOrCreateFieldCoordinates(fm, name="split coordinates")
        splitNodetemplate1 = nodes.createNodetemplate()
        splitNodetemplate1.defineField(splitCoordinates)
        splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        if useCubicHermiteThroughWall:
            splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
            if useCrossDerivatives:
                splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
                splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
                splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        splitNodetemplate2 = nodes.createNodetemplate()
        splitNodetemplate2.defineField(splitCoordinates)
        splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_VALUE, 2)
        splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS1, 2)
        splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS2, 2)
        if useCrossDerivatives:
            splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 2)
        if useCubicHermiteThroughWall:
            splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS3, 2)
            if useCrossDerivatives:
                splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 2)
                splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 2)
                splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 2)

        nodeIter = nodes.createNodeiterator()
        node = nodeIter.next()
        while node.isValid():
            if not markerPoints.containsNode(node):
                cache.setNode(node)
                identifier = node.getIdentifier()
                marginNode = identifier in nodesOnSplitMargin
                x = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)[1]
                d1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)[1]
                d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)[1]
                if useCrossDerivatives:
                    d1d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, 3)[1]
                if useCubicHermiteThroughWall:
                    d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)[1]
                    if useCrossDerivatives:
                        d1d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, 3)[1]
                        d2d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, 3)[1]
                        d1d2d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, 3)[1]

                node.merge(splitNodetemplate2 if marginNode else splitNodetemplate1)
                versionCount = 2 if marginNode else 1
                for vn in range(versionCount):
                    splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, vn + 1, x)
                    splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, vn + 1, d1)
                    splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, vn + 1, d2)
                    if useCrossDerivatives:
                        splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, vn + 1, d1d2)
                    if useCubicHermiteThroughWall:
                        splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, vn + 1, d3)
                        if useCrossDerivatives:
                            splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, vn + 1, d1d3)
                            splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, vn + 1, d2d3)
                            splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, vn + 1, d1d2d3)

            node = nodeIter.next()

        elementIter = mesh.createElementiterator()
        element = elementIter.next()
        splitElementtemplate1 = mesh.createElementtemplate()
        splitElementtemplate2 = mesh.createElementtemplate()

        count = 0
        elementsInOstium = elementsCountAroundEso * elementsAlongEsophagus * elementsThroughEsophagusWall
        closedLoopElementId = nextElementIdentifier - elementsCountAroundEso * elementsCountAcrossCardia - \
                              elementsCountAroundDuod * elementsCountThroughWall * (elementsAlongCardiaToDuod + 1)

        allValueLabels = [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                           Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3,
                           Node.VALUE_LABEL_D2_DS2DS3, Node.VALUE_LABEL_D3_DS1DS2DS3]

        while element.isValid():
            eft = element.getElementfieldtemplate(coordinates, -1)
            nodeIdentifiers = get_element_node_identifiers(element, eft)
            elementId = element.getIdentifier()
            marginDorsal = False
            for n in range(len(nodeIdentifiers)):
                marginElement = nodeIdentifiers[n] in nodesOnSplitMargin
                if marginElement:
                    count += 1
                    if count < 3 and (elementId <= elementsInOstium or elementId > closedLoopElementId):
                        marginDorsal = True
                    elif count >= 3 and (elementId <= elementsInOstium or elementId > closedLoopElementId):
                        if count == 4:
                            count = 0
                    elif elementsInOstium < elementId < elementsInOstium + len(xOuter[0]) + 1:
                        marginDorsal = True
                        count = 0
                    elif elementsInOstium + len(xOuter[0]) < elementId < elementsInOstium + len(xOuter[0]) * 2 + 1:
                        count = 0
                    elif count < 2 and elementId > elementsInOstium + 2 * (len(xOuter[0])):
                        marginDorsal = True
                    elif count >= 2 and elementId > elementsInOstium + 2 * (len(xOuter[0])):
                        if count == 2:
                            count = 0
                    break

            if marginDorsal:
                # Find nodes on margin to remap with version 2
                lnRemapV2 = []
                for n in range(len(nodeIdentifiers)):
                    if nodeIdentifiers[n] in nodesOnSplitMargin:
                        lnRemapV2.append(n + 1)
                eft2 = eft
                remapEftNodeValueLabelsVersion(eft2, lnRemapV2, allValueLabels, 2)

                result1 = splitElementtemplate2.defineField(splitCoordinates, -1, eft2)
                result2 = element.merge(splitElementtemplate2)
                element.setNodesByIdentifier(eft2, nodeIdentifiers)
                if eft2.getNumberOfLocalScaleFactors() > 0:
                    element.setScaleFactor(eft2, 1, -1.0)
            else:
                result1 = splitElementtemplate1.defineField(splitCoordinates, -1, eft)
                result2 = element.merge(splitElementtemplate1)
                element.setNodesByIdentifier(eft, nodeIdentifiers)
                if eft.getNumberOfLocalScaleFactors() == 1:
                    element.setScaleFactors(eft, [-1.0])

            element = elementIter.next()

        fm.endChange()

        return allAnnotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        refineElementsCountAround = options['Refine number of elements surface']
        refineElementsCountAlong = options['Refine number of elements surface']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong,
                                                       refineElementsCountThroughWall)
        return

    @classmethod
    def defineFaceAnnotations(cls, region, options, annotationGroups):
        '''
        Add face annotation groups from the highest dimension mesh.
        Must have defined faces and added subelements for highest dimension groups.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for top-level elements.
        New face annotation groups are appended to this list.
        '''

        limitingRidge = options['Limiting ridge']

        if limitingRidge:
            fm = region.getFieldmodule()
            fundusGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("fundus of stomach"))
            bodyGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("body of stomach"))
            cardiaGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("cardia of stomach"))
            antrumGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("pyloric antrum"))

            limitingRidgeGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term("forestomach-glandular stomach junction"))
            innerLimitingRidgeGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term("forestomach-glandular stomach junction on inner wall"))
            outerLimitingRidgeGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term("forestomach-glandular stomach junction on outer wall"))

            mesh2d = fm.findMeshByDimension(2)
            is_fundus = fundusGroup.getGroup()
            is_body = bodyGroup.getGroup()
            is_cardia = cardiaGroup.getGroup()
            is_antrum = antrumGroup.getGroup()
            is_limitingRidgeBody = fm.createFieldAnd(is_fundus, is_body)
            is_limitingRidgeCardia = fm.createFieldAnd(is_body, is_cardia)
            is_limitingRidgeAntrum = fm.createFieldAnd(is_antrum, is_cardia)
            is_limitingRidgeBodyCardia = fm.createFieldOr(is_limitingRidgeBody, is_limitingRidgeCardia)
            is_limitingRidgeAntrumCardia = fm.createFieldOr(is_limitingRidgeAntrum, is_limitingRidgeCardia)
            is_limitingRidge = fm.createFieldOr(is_limitingRidgeBodyCardia, is_limitingRidgeAntrumCardia)
            limitingRidgeGroup.getMeshGroup(mesh2d).addElementsConditional(is_limitingRidge)

            mesh1d = fm.findMeshByDimension(1)
            is_exterior = fm.createFieldIsExterior()
            is_exterior_face_outer = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
            is_exterior_face_inner = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))

            is_limitingRidgeInner = fm.createFieldAnd(is_limitingRidge, is_exterior_face_inner)
            innerLimitingRidgeGroup.getMeshGroup(mesh1d).addElementsConditional(is_limitingRidgeInner)

            is_limitingRidgeOuter = fm.createFieldAnd(is_limitingRidge, is_exterior_face_outer)
            outerLimitingRidgeGroup.getMeshGroup(mesh1d).addElementsConditional(is_limitingRidgeOuter)

        serosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                         get_stomach_term("serosa of stomach"))
        mucosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                         get_stomach_term("mucosa of stomach"))

        outerSplitMarginGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                         get_stomach_term("serosa split margin"))
        esoSplitMarginGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                              get_stomach_term("esophagus on serosa split margin"))
        fundusSplitMarginGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                              get_stomach_term("fundus on serosa split margin"))
        bodySplitMarginGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                              get_stomach_term("body on serosa split margin"))
        antrumSplitMarginGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                     get_stomach_term("antrum on serosa split margin"))
        pylorusSplitMarginGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                              get_stomach_term("pylorus on serosa split margin"))

        serosaDorsalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                         get_stomach_term("serosa of dorsal stomach"))
        serosaVentralGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                         get_stomach_term("serosa of ventral stomach"))

        stomachGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("stomach"))
        bodyGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("body of stomach"))
        esoGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("esophagus"))
        fundusGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("fundus of stomach"))
        antrumGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("pyloric antrum"))
        pylorusGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("pylorus"))
        esoGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("esophagus"))
        esoGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("esophagus"))
        dorsalGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("dorsal stomach"))
        ventralGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("ventral stomach"))

        fm = region.getFieldmodule()
        mesh2d = fm.findMeshByDimension(2)
        is_stomach = stomachGroup.getGroup()
        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_outer = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_exterior_face_inner = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))
        is_serosa = fm.createFieldAnd(is_stomach, is_exterior_face_outer)
        is_mucosa = fm.createFieldAnd(is_stomach, is_exterior_face_inner)
        serosaGroup.getMeshGroup(mesh2d).addElementsConditional(is_serosa)
        mucosaGroup.getMeshGroup(mesh2d).addElementsConditional(is_mucosa)
        is_dorsal = dorsalGroup.getGroup()
        is_ventral = ventralGroup.getGroup()
        is_dorsalSerosa = fm.createFieldAnd(is_dorsal, is_serosa)
        is_ventralSerosa = fm.createFieldAnd(is_ventral, is_serosa)
        serosaDorsalGroup.getMeshGroup(mesh2d).addElementsConditional(is_dorsalSerosa)
        serosaVentralGroup.getMeshGroup(mesh2d).addElementsConditional(is_ventralSerosa)

        is_margin = fm.createFieldAnd(is_dorsal, is_ventral)
        mesh1d = fm.findMeshByDimension(1)
        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_outer = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_marginOuter = fm.createFieldAnd(is_margin, is_exterior_face_outer)
        outerSplitMarginGroup.getMeshGroup(mesh1d).addElementsConditional(is_marginOuter)

        is_eso = esoGroup.getGroup()
        is_esoMargin = fm.createFieldAnd(is_marginOuter, is_eso)
        esoSplitMarginGroup.getMeshGroup(mesh1d).addElementsConditional(is_esoMargin)

        is_fundus = fundusGroup.getGroup()
        is_fundusMargin = fm.createFieldAnd(is_marginOuter, is_fundus)
        fundusSplitMarginGroup.getMeshGroup(mesh1d).addElementsConditional(is_fundusMargin)

        is_body = bodyGroup.getGroup()
        is_bodyMargin = fm.createFieldAnd(is_marginOuter, is_body)
        bodySplitMarginGroup.getMeshGroup(mesh1d).addElementsConditional(is_bodyMargin)

        is_antrum = antrumGroup.getGroup()
        is_antrumMargin = fm.createFieldAnd(is_marginOuter, is_antrum)
        antrumSplitMarginGroup.getMeshGroup(mesh1d).addElementsConditional(is_antrumMargin)

        is_pylorus = pylorusGroup.getGroup()
        is_pylorusMargin = fm.createFieldAnd(is_marginOuter, is_pylorus)
        pylorusSplitMarginGroup.getMeshGroup(mesh1d).addElementsConditional(is_pylorusMargin)


def findClosestPositionAndDerivativeOnTrackSurface(x, nx, trackSurface, nxProportion1, elementsCountAlongTrackSurface):
    """
    Find the closest position and derivative around the tracksurface of a point sitting near the fundus of stomach.
    Use a startPosition to improve the search for nearest position on the track surface as the fundus has a curved
    and complex track surface.
    :param x: coordinates of point of interest
    :param nx: coordinates of points along curve where point of interest lies.
    :param trackSurface: track surface where point sits
    :param nxProportion1: proportion around track surface of curve
    :param elementsCountAlongTrackSurface: number of elements along track surface
    :return: position and derivative of point around track surface
    """
    closestIdxOnNx = interp.getNearestPointIndex(nx, x)
    closestPositionToPoint = trackSurface.createPositionProportion(nxProportion1, closestIdxOnNx / elementsCountAlongTrackSurface)
    xPosition = trackSurface.findNearestPosition(x, closestPositionToPoint)
    d = trackSurface.evaluateCoordinates(xPosition, derivatives=True)[1]

    return xPosition, d

def getSmoothedSampledPointsOnTrackSurface(trackSurface, startProportion1, startProportion2, endProportion1,
                                           endProportion2, elementsOut, startDerivative = None, endDerivative = None,
                                           startDerivativeMagnitude = None, endDerivativeMagnitude = None,
                                           curveMode = 1):
    """
    Create smoothly spaced out hermite curve points between two points a and b on the surface,
    each defined by their proportions over the surface in directions 1 and 2.
    :param trackSurface: track surface
    :param startProportion1, startProportion2: proportion of start point in direction around and along track surface
    :param endProportion1, endProportion2: proportion of end point in direction around and along track surface
    :param elementsOut: number of elements out
    :param startDerivative, endDerivative: optional derivative vectors in 3-D world coordinates
        to match at the start and end of the curves. If omitted, fits in with other derivative or is
        in a straight line from a to b
    :param derivativeMagnitudeStart, derivativeMagnitudeEnd: optional magnitude of derivatives to match at the start and
        end of the curves
    :return: coordinates and derivative of sampled points
    """

    mx, md2, md1, md3, mProportions = \
        trackSurface.createHermiteCurvePoints(startProportion1, startProportion2, endProportion1, endProportion2,
                                              elementsOut, startDerivative, endDerivative, curveMode)

    xSampled, dSampled = trackSurface.resampleHermiteCurvePointsSmooth(mx, md2, md1, md3, mProportions,
                                                                       startDerivativeMagnitude,
                                                                       endDerivativeMagnitude)[0:2]
    return xSampled, dSampled

def smoothD1Around(xAround, d1Around):
    """
    Rearrange points around so that the group of points starts from left side of the annulus to the greater curvature
    and ends on the right side of the annulus. This put points in consecutive order for derivative smoothing.
    The smoothed derivatives are then re-arranged such that it starts from the greater curvature and goes to the right
    side of the annulus, followed by the left side of the annulus and closing the loop back at the greater curvature.
    :param xAround: points around a loop joining to the annulus.
    :param d1Around: derivative of points.
    :return: smoothed derivatives in the original order.
    """
    xLoop = xAround[int(len(xAround) * 0.5 + 1): ] + xAround[: int(len(xAround) * 0.5 + 1)]
    d1Loop = d1Around[int(len(d1Around) * 0.5 + 1): ] + d1Around[: int(len(d1Around) * 0.5 + 1)]
    d1LoopSmooth = interp.smoothCubicHermiteDerivativesLine(xLoop, d1Loop, fixStartDerivative=True,
                                                                            fixEndDerivative=True)
    # Rearrange to correct order
    d1Around = d1LoopSmooth[int(len(xAround) * 0.5) : ] + d1LoopSmooth[: int(len(xAround) * 0.5) : ]

    return d1Around

def findD1CurvatureAround(xAround, d1Around, normsAround):
    """
    Rearrange points around so that the group of points starts from left side of the annulus to the greater curvature
    and ends on the right side of the annulus. This put points in consecutive order for calculating curvature.
    The calculated curvatures are then re-arranged such that it starts from the greater curvature and goes to the right
    side of the annulus, followed by the left side of the annulus and closing the loop back at the greater curvature.
    :param xAround: points around a loop joining to the annulus.
    :param d1Around: derivative of points.
    :param normsAround: radial normal at each point.
    :return: curvature in the original order.
    """
    xLoop = xAround[int(len(xAround) * 0.5 + 1):] + xAround[: int(len(xAround) * 0.5 + 1)]
    d1Loop = d1Around[int(len(d1Around) * 0.5 + 1):] + d1Around[: int(len(d1Around) * 0.5 + 1)]
    normsLoop = normsAround[int(len(d1Around) * 0.5 + 1):] + normsAround[: int(len(d1Around) * 0.5 + 1)]
    curvature = findCurvatureAlongLine(xLoop, d1Loop, normsLoop)
    # Rearrange to correct order
    d1CurvatureAround = curvature[int(len(xAround) * 0.5):] + curvature[: int(len(xAround) * 0.5):]

    return d1CurvatureAround

def findDerivativeBetweenPoints(v1, v2):
    """
    Find vector difference between two points and rescale vector difference using cubic hermite arclength
    between the points to derive the derivative between the points.
    :param v1: start vector
    :param v2: end vector
    :return: derivative of between v1 and v2
    """
    d = [v2[c] - v1[c] for c in range(3)]
    arcLengthAround = interp.computeCubicHermiteArcLength(v1, d, v2, d, True)
    d = [c * arcLengthAround for c in vector.normalise(d)]

    return d

def findCurvatureAroundLoop(nx, nd, radialVectors):
    """
    Calculate curvature for points lying along a loop.
    :param nx: points on loop
    :param nd: derivative of points on loop
    :param radialVectors: radial direction, assumed normal to curve tangent at point
    :return: curvatures along points on loop
    """
    curvature = []
    for n in range(len(nx)):
        prevIdx = n - 1 if n > 0 else -1
        nextIdx = n + 1 if n < len(nx) - 1 else 0
        kappam = interp.getCubicHermiteCurvature(nx[prevIdx], nd[prevIdx], nx[n], nd[n], radialVectors[n], 1.0)
        kappap = interp.getCubicHermiteCurvature(nx[n], nd[n], nx[nextIdx], nd[nextIdx], radialVectors[n], 0.0)
        curvature.append(0.5 * (kappam + kappap))

    return curvature

def findCurvatureAlongLine(nx, nd, radialVectors):
    """
    Calculate curvature for points lying along a line.
    :param nx: points on line
    :param nd: derivative of points on line
    :param radialVectors: radial direction, assumed normal to curve tangent at point
    :return: curvatures along points on line
    """
    curvature = []
    for n in range(len(nx)):
        if n == 0:
            curvature.append(interp.getCubicHermiteCurvature(nx[n], nd[n], nx[n + 1], nd[n + 1], radialVectors[n], 0.0))
        elif n == len(nx) - 1:
            curvature.append(interp.getCubicHermiteCurvature(nx[n - 1], nd[n - 1], nx[n], nd[n], radialVectors[n], 1.0))
        else:
            curvature.append(0.5 * (
                        interp.getCubicHermiteCurvature(nx[n], nd[n], nx[n + 1], nd[n + 1], radialVectors[n], 0.0) +
                        interp.getCubicHermiteCurvature(nx[n - 1], nd[n - 1], nx[n], nd[n], radialVectors[n], 1.0)))

    return curvature
