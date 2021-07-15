"""
Generates body coordinates using a solid cylinder of all cube elements,
 with variable numbers of elements in major, minor, shell and for each section of abdomen, thorax, neck and head.
"""

from __future__ import division
import math
import copy
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup,\
    findOrCreateFieldStoredString, findOrCreateFieldStoredMeshLocation, findOrCreateFieldNodeGroup
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, CylinderEnds, CylinderCentralPath, createEllipsePerimeter
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from opencmiss.zinc.node import Node
from opencmiss.zinc.element import Element
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm, mergeAnnotationGroups
from scaffoldmaker.annotation.body_terms import get_body_term
from scaffoldmaker.annotation import heart_terms, bladder_terms, lung_terms, stomach_terms, brainstem_terms
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.field import Field, FieldFindMeshLocation
from opencmiss.utils.zinc.finiteelement import get_element_node_identifiers
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabelsVersion, remapEftNodeValueLabel, setEftScaleFactorIds
from scaffoldmaker.utils.vector import setMagnitude, crossproduct3
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite


class MeshType_3d_wholebody1(Scaffold_base):
    """
Generates body coordinates using a solid cylinder of all cube elements,
 with variable numbers of elements in major, minor, shell and for each section of abdomen, thorax, neck and head.
    """
    cylinder1Settings = {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.5,
                'Number of elements': 1
            }

    axis1 = [0, 0, 1]
    axis2 = [1, 0, 0]
    axis3 = [0, 1, 0]
    centralPathDefaultScaffoldPackages = {
        'Default': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': cylinder1Settings,
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[0.0, 0.0, 0.0], setMagnitude(axis1, cylinder1Settings['Length']), setMagnitude(axis2, 0.5), [0.0, 0.0, 0.0], setMagnitude(axis3, 0.5), [0.0, 0.0, 0.0]],
                    [setMagnitude(axis1, cylinder1Settings['Length']), setMagnitude(axis1, cylinder1Settings['Length']), setMagnitude(axis2, 0.5), [0.0, 0.0, 0.0], setMagnitude(axis3, 0.5), [0.0, 0.0, 0.0]]
                ])
        })
    }

    @staticmethod
    def getName():
        return '3D Whole Body 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human Coarse',
            'Human Fine',
            'Rat Coarse',
            'Rat Fine'
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if parameterSetName == 'Default':
            parameterSetName = 'Human Coarse'
        centralPathOption = cls.centralPathDefaultScaffoldPackages['Default']
        options = {}
        options['Base parameter set'] = parameterSetName
        options['Central path'] = copy.deepcopy(centralPathOption)
        options['Number of elements across major'] = 6
        options['Number of elements across minor'] = 6
        options['Number of elements across shell'] = 1
        options['Number of elements across transition'] = 1
        options['Number of elements in abdomen'] = 5
        options['Number of elements in thorax'] = 3
        options['Number of elements in neck'] = 1
        options['Number of elements in head'] = 2
        options['Shell thickness proportion'] = 0.2
        options['Discontinuity on the core boundary'] = True
        options['Lower half'] = False
        options['Use cross derivatives'] = False
        options['Refine'] = False
        options['Refine number of elements across major'] = 1
        options['Refine number of elements along'] = 1

        if 'Coarse' in parameterSetName:
            pass
        if 'Fine' in parameterSetName:
            options['Number of elements across major'] = 10
            options['Number of elements across minor'] = 10
            options['Number of elements across shell'] = 1
            options['Number of elements across transition'] = 1
            options['Number of elements in abdomen'] = 10
            options['Number of elements in thorax'] = 6
            options['Number of elements in neck'] = 2
            options['Number of elements in head'] = 4
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements across major',
            'Number of elements across minor',
            'Number of elements across shell',
            'Number of elements across transition',
            'Number of elements in abdomen',
            'Number of elements in thorax',
            'Number of elements in neck',
            'Number of elements in head',
            'Shell thickness proportion',
            'Discontinuity on the core boundary',
            'Refine',
            'Refine number of elements across major',
            'Refine number of elements along'
        ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [MeshType_1d_path1]
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
        dependentChanges = False

        if options['Number of elements across major'] < 4:
            options['Number of elements across major'] = 4
        if options['Number of elements across major'] % 2:
            options['Number of elements across major'] += 1

        if options['Number of elements across minor'] != options['Number of elements across major']:
            options['Number of elements across minor'] = options['Number of elements across major']
            dependentChanges = True

        if options['Number of elements across transition'] < 1:
            options['Number of elements across transition'] = 1

        Rcrit = min(options['Number of elements across major']-4, options['Number of elements across minor']-4)//2
        if options['Number of elements across shell'] + options['Number of elements across transition'] - 1 > Rcrit:
            dependentChanges = True
            options['Number of elements across shell'] = Rcrit
            options['Number of elements across transition'] = 1

        if options['Shell thickness proportion'] < 0.07 or options['Shell thickness proportion'] > 0.7:
            options['Shell thickness proportion'] = 2*options['Number of elements across shell']/options['Number of elements across major']

        if options['Number of elements in abdomen'] < 1:
            options['Number of elements in abdomen'] = 1
        if options['Number of elements in head'] < 1:
            options['Number of elements in head'] = 1
        if options['Number of elements in neck'] < 1:
            options['Number of elements in neck'] = 1
        if options['Number of elements in thorax'] < 1:
            options['Number of elements in thorax'] = 1

        return dependentChanges

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: List of AnnotationGroup
        """

        baseParameterSetName = options['Base parameter set']
        isHuman = 'Human' in baseParameterSetName
        isRat = 'Rat' in baseParameterSetName

        centralPath = options['Central path']
        full = not options['Lower half']
        elementsCountAcrossMajor = options['Number of elements across major']
        if not full:
            elementsCountAcrossMajor //= 2
        elementsCountAcrossMinor = options['Number of elements across minor']
        elementsCountAcrossShell = options['Number of elements across shell']
        elementsCountAcrossTransition = options['Number of elements across transition']
        elementsCountAlongAbdomen = options['Number of elements in abdomen']
        elementsCountAlongHead = options['Number of elements in head']
        elementsCountAlongNeck = options['Number of elements in neck']
        elementsCountAlongThorax = options['Number of elements in thorax']
        shellRadiusProportion = options['Shell thickness proportion']
        shellProportion = 1/(1/shellRadiusProportion-1)*(elementsCountAcrossMajor/2/elementsCountAcrossShell - 1)
        discontinuity = options['Discontinuity on the core boundary']
        useCrossDerivatives = options['Use cross derivatives']

        elementsCountAlong = elementsCountAlongAbdomen + elementsCountAlongThorax + elementsCountAlongNeck + elementsCountAlongHead

        fieldmodule = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fieldmodule)
        mesh = fieldmodule.findMeshByDimension(3)

        bodyGroup = AnnotationGroup(region, get_body_term("body"))
        coreGroup = AnnotationGroup(region, get_body_term("core"))
        non_coreGroup = AnnotationGroup(region, get_body_term("non core"))
        abdomenGroup = AnnotationGroup(region, get_body_term("abdomen"))
        thoraxGroup = AnnotationGroup(region, get_body_term("thorax"))
        neckGroup = AnnotationGroup(region, get_body_term("neck core"))
        headGroup = AnnotationGroup(region, get_body_term("head core"))
        annotationGroups = [bodyGroup, coreGroup, non_coreGroup, abdomenGroup, thoraxGroup, neckGroup, headGroup]

        cylinderCentralPath = CylinderCentralPath(region, centralPath, elementsCountAlong)

        cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL

        base = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossShell,
                            elementsCountAcrossTransition, shellProportion,
                            [0.0, 0.0, 0.0], cylinderCentralPath.alongAxis[0], cylinderCentralPath.majorAxis[0],
                            cylinderCentralPath.minorRadii[0])
        cylinder1 = CylinderMesh(fieldmodule, coordinates, elementsCountAlong, base,
                                 cylinderShape=cylinderShape,
                                 cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False)

        # body coordinates
        bodyCoordinates = findOrCreateFieldCoordinates(fieldmodule, name="body coordinates")
        tmp_region = region.createRegion()
        tmp_fieldmodule = tmp_region.getFieldmodule()
        tmp_body_coordinates = findOrCreateFieldCoordinates(tmp_fieldmodule, name="body coordinates")
        tmp_cylinder = CylinderMesh(tmp_fieldmodule, tmp_body_coordinates, elementsCountAlong, base,
                                 cylinderShape=cylinderShape,
                                 cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False)
        sir = tmp_region.createStreaminformationRegion()
        srm = sir.createStreamresourceMemory()
        tmp_region.write(sir)
        result, buffer = srm.getBuffer()
        sir = region.createStreaminformationRegion()
        srm = sir.createStreamresourceMemoryBuffer(buffer)
        region.read(sir)

        del srm
        del sir
        del tmp_body_coordinates
        del tmp_fieldmodule
        del tmp_region

        # Groups of different parts of the body
        is_body = fieldmodule.createFieldConstant(1)
        bodyMeshGroup = bodyGroup.getMeshGroup(mesh)
        bodyMeshGroup.addElementsConditional(is_body)

        coreMeshGroup = coreGroup.getMeshGroup(mesh)

        # core group
        e1a = elementsCountAcrossShell
        e1z = elementsCountAcrossMinor - elementsCountAcrossShell - 1
        e2a = elementsCountAcrossShell
        e2b = e2a + elementsCountAcrossTransition
        e2z = elementsCountAcrossMajor - elementsCountAcrossShell - 1
        e2y = e2z - elementsCountAcrossTransition
        e1oc = elementsCountAcrossMinor - 2*elementsCountAcrossShell - 2*elementsCountAcrossTransition
        e2oc = elementsCountAcrossMajor - 2*elementsCountAcrossShell - 2*elementsCountAcrossTransition
        e2oCore = e2oc * e1oc + 2 * elementsCountAcrossTransition * (e2oc + e1oc)
        elementsCountAround = cylinder1.getElementsCountAround()
        e2oShell = elementsCountAround * elementsCountAcrossShell
        e2o = e2oCore + e2oShell
        elementId = cylinder1.getElementIdentifiers()
        for e3 in range(elementsCountAlong):
            for e2 in range(elementsCountAcrossMajor):
                for e1 in range(elementsCountAcrossMinor):
                    coreElement = ((e2 >= e2a) and (e2 <= e2z)) and ((e1 >= e1a) and (e1 <= e1z))
                    if coreElement:
                        elementIdentifier = elementId[e3][e2][e1]
                        if elementIdentifier:
                            element = mesh.findElementByIdentifier(elementIdentifier)
                            coreMeshGroup.addElement(element)

        is_non_core = fieldmodule.createFieldNot(coreGroup.getGroup())
        non_coreMeshGroup = non_coreGroup.getMeshGroup(mesh)
        non_coreMeshGroup.addElementsConditional(is_non_core)

        abdomenMeshGroup = abdomenGroup.getMeshGroup(mesh)
        thoraxMeshGroup = thoraxGroup.getMeshGroup(mesh)
        neckMeshGroup = neckGroup.getMeshGroup(mesh)
        headMeshGroup = headGroup.getMeshGroup(mesh)
        meshGroups = [abdomenMeshGroup, thoraxMeshGroup, neckMeshGroup, headMeshGroup]

        abdomenRange = [1, elementsCountAlongAbdomen*e2o]
        thoraxRange = [abdomenRange[1]+1, abdomenRange[1]+elementsCountAlongThorax*e2o]
        neckRange = [thoraxRange[1]+1, thoraxRange[1] + elementsCountAlongNeck*e2o]
        headRange = [neckRange[1]+1, elementsCountAlong*e2o]
        groupsRanges = [abdomenRange, thoraxRange, neckRange, headRange]

        totalElements = e2o*elementsCountAlong
        for elementIdentifier in range(1, totalElements+1):
            element = mesh.findElementByIdentifier(elementIdentifier)
            if coreMeshGroup.containsElement(element):
                ri = 0
                for groupRange in groupsRanges:
                    if (elementIdentifier >= groupRange[0]) and (elementIdentifier <= groupRange[1]):
                        meshGroups[ri].addElement(element)
                        break
                    ri += 1

        if discontinuity:
            # create discontinuity in d3 on the core boundary
            nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
            elementtemplate = mesh.createElementtemplate()
            undefineNodetemplate = nodes.createNodetemplate()
            undefineNodetemplate.undefineField(coordinates)
            nodetemplate = nodes.createNodetemplate()
            fieldcache = fieldmodule.createFieldcache()
            with ChangeManager(fieldmodule):
                localNodeIndexes = [1, 2, 3, 4]
                valueLabel = Node.VALUE_LABEL_D_DS3
                for e3 in range(elementsCountAlong):
                    for e2 in range(elementsCountAcrossMajor):
                        for e1 in range(elementsCountAcrossMinor):
                            regularRowElement = (((e2 >= e2b) and (e2 <= e2y)) and ((e1 == e1a - 1) or (e1 == e1z + 1)))
                            non_coreFirstLayerElement = (e2 == e2a - 1) or regularRowElement or (e2 == e2z + 1)
                            elementIdentifier = elementId[e3][e2][e1]
                            if elementIdentifier and non_coreFirstLayerElement:
                                element = mesh.findElementByIdentifier(elementIdentifier)
                                eft = element.getElementfieldtemplate(coordinates, -1)
                                nodeIds = get_element_node_identifiers(element, eft)
                                for localNodeIndex in localNodeIndexes:
                                    node = element.getNode(eft, localNodeIndex)
                                    nodetemplate.defineFieldFromNode(coordinates, node)
                                    versionsCount = nodetemplate.getValueNumberOfVersions(coordinates, -1, valueLabel)
                                    if versionsCount == 1:
                                        fieldcache.setNode(node)
                                        result0, x = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                                        result0, d1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                                        result0, d2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                                        result0, d3 = coordinates.getNodeParameters(fieldcache, -1, valueLabel, 1, 3)
                                        result1 = node.merge(undefineNodetemplate)
                                        result2 = nodetemplate.setValueNumberOfVersions(coordinates, -1, valueLabel, 2)
                                        result3 = node.merge(nodetemplate)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1, valueLabel, 1, d3)
                                        result5 = coordinates.setNodeParameters(fieldcache, -1, valueLabel, 2, d3)
                                remapEftNodeValueLabelsVersion(eft, localNodeIndexes, [valueLabel], 2)
                                result1 = elementtemplate.defineField(coordinates, -1, eft)
                                result2 = element.merge(elementtemplate)
                                result3 = element.setNodesByIdentifier(eft, nodeIds)
        else:
            fieldcache = fieldmodule.createFieldcache()

        # Annotation fiducial point
        markerGroup = findOrCreateFieldGroup(fieldmodule, "marker")
        markerName = findOrCreateFieldStoredString(fieldmodule, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fieldmodule, mesh, name="marker_location")
        markerBodyCoordinates = findOrCreateFieldCoordinates(fieldmodule, name="marker_body_coordinates")
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)
        markerTemplateInternal.defineField(markerBodyCoordinates)
        #
        middleLeft = elementsCountAcrossMinor//2
        topElem = elementsCountAcrossMajor - 1
        middleRight = middleLeft - 1
        neckFirstElem = elementsCountAlongAbdomen+elementsCountAlongThorax
        thoraxFirstElem = elementsCountAlongAbdomen
        middleDown = elementsCountAcrossMajor//2 - 1

        # organ landmarks groups
        apexOfHeart = heart_terms.get_heart_term('apex of heart')
        leftAtriumEpicardiumVenousMidpoint = heart_terms.get_heart_term('left atrium epicardium venous midpoint')
        rightAtriumEpicardiumVenousMidpoint = heart_terms.get_heart_term('right atrium epicardium venous midpoint')
        apexOfUrinaryBladder = bladder_terms.get_bladder_term('apex of urinary bladder')
        leftUreterJunctionWithBladder = bladder_terms.get_bladder_term('left ureter junction with bladder')
        rightUreterJunctionWithBladder = bladder_terms.get_bladder_term('right ureter junction with bladder')
        urethraJunctionWithBladderDorsal = bladder_terms.get_bladder_term('urethra junction of dorsal bladder neck')
        urethraJunctionWithBladderVentral = bladder_terms.get_bladder_term('urethra junction of ventral bladder neck')
        gastroesophagalJunctionOnLesserCurvature = stomach_terms.get_stomach_term('gastro-esophagal junction on lesser curvature')
        limitingRidgeOnGreaterCurvature = stomach_terms.get_stomach_term('limiting ridge on greater curvature')
        pylorusOnGreaterCurvature = stomach_terms.get_stomach_term('pylorus on greater curvature')
        duodenumOnGreaterCurvature = stomach_terms.get_stomach_term('duodenum on greater curvature')
        junctionBetweenFundusAndBodyOnGreaterCurvature = stomach_terms.get_stomach_term("junction between fundus and body on greater curvature")
        apexOfLeftLung = lung_terms.get_lung_term('apex of left lung')
        ventralBaseOfLeftLung = lung_terms.get_lung_term('ventral base of left lung')
        dorsalBaseOfLeftLung = lung_terms.get_lung_term('dorsal base of left lung')
        apexOfRightLung = lung_terms.get_lung_term('apex of right lung')
        ventralBaseOfRightLung = lung_terms.get_lung_term('ventral base of right lung')
        dorsalBaseOfRightLung = lung_terms.get_lung_term('dorsal base of right lung')
        laterodorsalTipOfMiddleLobeOfRightLung = lung_terms.get_lung_term('laterodorsal tip of middle lobe of right lung')
        apexOfRightLungAccessoryLobe = lung_terms.get_lung_term('apex of right lung accessory lobe')
        ventralBaseOfRightLungAccessoryLobe = lung_terms.get_lung_term('ventral base of right lung accessory lobe')
        dorsalBaseOfRightLungAccessoryLobe = lung_terms.get_lung_term('dorsal base of right lung accessory lobe')
        medialBaseOfLeftLung = lung_terms.get_lung_term("medial base of left lung")
        medialBaseOfRightLung = lung_terms.get_lung_term("medial base of right lung")
        brainstemDorsalMidlineCaudalPoint = brainstem_terms.get_brainstem_term('brainstem dorsal midline caudal point')
        brainstemDorsalMidlineCranialPoint = brainstem_terms.get_brainstem_term('brainstem dorsal midline cranial point')
        brainstemVentralMidlineCaudalPoint = brainstem_terms.get_brainstem_term('brainstem ventral midline caudal point')
        brainstemVentralMidlineCranialPoint = brainstem_terms.get_brainstem_term('brainstem ventral midline cranial point')

        # marker coordinates. In future we want to have only one table for all species.
        if isRat:
            bodyMarkerPoints = [
                {"group": ("left hip joint", ''), "x": [0.367, 0.266, 0.477]},
                {"group": ("right hip joint", ''), "x": [-0.367, 0.266, 0.477]},
                {"group": ("left shoulder joint", ''), "x": [0.456, -0.071, 2.705]},
                {"group": ("right shoulder joint", ''), "x": [-0.456, -0.071, 2.705]},
                {"group": ("along left femur", ''), "x": [0.456, 0.07, 0.633]},
                {"group": ("along right femur", ''), "x": [-0.456, 0.07, 0.633]},
                {"group": ("along left humerus", ''), "x": [0.423, -0.173, 2.545]},
                {"group": ("along right humerus", ''), "x": [-0.423, -0.173, 2.545]},
                {"group": apexOfUrinaryBladder, "x": [-0.124, -0.383, 0.434]},
                {"group": leftUreterJunctionWithBladder, "x": [-0.111, -0.172, 0.354]},
                {"group": rightUreterJunctionWithBladder, "x": [-0.03, -0.196, 0.363]},
                {"group": urethraJunctionWithBladderDorsal, "x": [-0.03, -0.26, 0.209]},
                {"group": urethraJunctionWithBladderVentral, "x": [-0.037, -0.304, 0.203]},
                {"group": brainstemDorsalMidlineCaudalPoint, "x": [-0.032, 0.418, 2.713]},
                {"group": brainstemDorsalMidlineCranialPoint, "x": [-0.017, 0.203, 2.941]},
                {"group": brainstemVentralMidlineCaudalPoint, "x": [-0.028, 0.388, 2.72]},
                {"group": brainstemVentralMidlineCranialPoint, "x": [-0.019, 0.167, 2.95]},
                {"group": apexOfHeart, "x": [0.096, -0.128, 1.601]},
                {"group": leftAtriumEpicardiumVenousMidpoint, "x": [0.127, -0.083, 2.079]},
                {"group": rightAtriumEpicardiumVenousMidpoint, "x": [0.039, -0.082, 2.075]},
                {"group": apexOfLeftLung, "x": [0.172, -0.175, 2.337]},
                {"group": ventralBaseOfLeftLung, "x": [0.274, -0.285, 1.602]},
                {"group": dorsalBaseOfLeftLung, "x": [0.037, 0.31, 1.649]},
                {"group": apexOfRightLung, "x": [-0.086, -0.096, 2.311]},
                {"group": ventralBaseOfRightLung, "x": [0.14, -0.357, 1.662]},
                {"group": dorsalBaseOfRightLung, "x": [-0.054, 0.304, 1.667]},
                {"group": laterodorsalTipOfMiddleLobeOfRightLung, "x": [-0.258, -0.173, 2.013]},
                {"group": apexOfRightLungAccessoryLobe, "x": [0.041, -0.063, 1.965]},
                {"group": ventralBaseOfRightLungAccessoryLobe, "x": [0.143, -0.356, 1.66]},
                {"group": dorsalBaseOfRightLungAccessoryLobe, "x": [0.121, -0.067, 1.627]},
                {"group": gastroesophagalJunctionOnLesserCurvature, "x": [0.12, 0.009, 1.446]},
                {"group": limitingRidgeOnGreaterCurvature, "x": [0.318, 0.097, 1.406]},
                {"group": pylorusOnGreaterCurvature, "x": [0.08, -0.111, 1.443]},
                {"group": duodenumOnGreaterCurvature, "x": [0.029, -0.138, 1.481]},
            ]
        elif isHuman:
            bodyMarkerPoints = [
                {"group": urethraJunctionWithBladderDorsal, "x": [-0.0071, -0.2439, 0.1798]},
                {"group": urethraJunctionWithBladderVentral, "x": [-0.007, -0.2528, 0.1732]},
                {"group": leftUreterJunctionWithBladder, "x": [0.1074, 0.045, 0.1728]},
                {"group": rightUreterJunctionWithBladder, "x": [-0.1058, 0.0533, 0.1701]},
                {"group": apexOfUrinaryBladder, "x": [0.005, 0.1286, 0.1264]},
                {"group": brainstemDorsalMidlineCaudalPoint, "x": [0.0068, 0.427, 2.7389]},
                {"group": brainstemDorsalMidlineCranialPoint, "x": [0.008, -0.0231, 3.0778]},
                {"group": brainstemVentralMidlineCaudalPoint, "x": [0.0054, 0.3041, 2.7374]},
                {"group": brainstemVentralMidlineCranialPoint, "x": [0.0025, -0.2308, 3.091]},
                {"group": apexOfHeart, "x": [0.1373, -0.1855, 1.421]},
                {"group": leftAtriumEpicardiumVenousMidpoint, "x": [0.0024, 0.1452, 1.8022]},
                {"group": rightAtriumEpicardiumVenousMidpoint, "x": [-0.0464, 0.0373, 1.7491]},
                {"group": apexOfLeftLung, "x": [0.0655, -0.0873, 2.3564]},
                {"group": apexOfRightLung, "x": [-0.088, -0.0363, 2.3518]},
                {"group": laterodorsalTipOfMiddleLobeOfRightLung, "x": [-0.2838, -0.0933, 1.9962]},
                {"group": ventralBaseOfLeftLung, "x": [0.219, -0.2866, 1.4602]},
                {"group": medialBaseOfLeftLung, "x": [0.0426, -0.0201, 1.4109]},
                {"group": ventralBaseOfRightLung, "x": [-0.2302, -0.2356, 1.3926]},
                {"group": medialBaseOfRightLung, "x": [-0.0363, 0.0589, 1.3984]},
                {"group": dorsalBaseOfLeftLung, "x": [0.1544, 0.2603, 1.3691]},
                {"group": dorsalBaseOfRightLung, "x": [0.0369, -0.2524, 0.912]},
                {"group": gastroesophagalJunctionOnLesserCurvature, "x": [-0.0062, -0.3259, 0.8586]},
                {"group": pylorusOnGreaterCurvature, "x": [-0.0761, -0.3189, 0.8663]},
                {"group": duodenumOnGreaterCurvature, "x": [-0.1599, 0.1601, 1.1939]},
                {"group": junctionBetweenFundusAndBodyOnGreaterCurvature, "x": [0.1884, -0.1839, 0.9639]},
            ]

        nodeIdentifier = cylinder1._endNodeIdentifier
        findMarkerLocation = fieldmodule.createFieldFindMeshLocation(markerBodyCoordinates, bodyCoordinates, mesh)
        findMarkerLocation.setSearchMode(FieldFindMeshLocation.SEARCH_MODE_EXACT)
        for bodyMarkerPoint in bodyMarkerPoints:
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            fieldcache.setNode(markerPoint)
            markerBodyCoordinates.assignReal(fieldcache, bodyMarkerPoint["x"])
            markerName.assignString(fieldcache, bodyMarkerPoint["group"][0])

            element, xi = findMarkerLocation.evaluateMeshLocation(fieldcache, 3)
            markerLocation.assignMeshLocation(fieldcache, element, xi)
            nodeIdentifier += 1





        # Join upper_limb
        upperlimb_along_proportion = 0.6
        upperlimb_elementCountAlong = round(upperlimb_along_proportion * cylinder1._elementsCountAlong)
        upperlimb_left_nid = cylinder1._shield.nodeId [upperlimb_elementCountAlong][0][cylinder1._elementsCountAcrossMinor//2]
        upperlimb_right_nid = cylinder1._shield.nodeId[upperlimb_elementCountAlong][-1][cylinder1._elementsCountAcrossMinor//2]

        circumferential_element_length = cylinder1._shield.px[upperlimb_elementCountAlong][0][cylinder1._elementsCountAcrossMinor//2][1] - cylinder1._shield.px[upperlimb_elementCountAlong][0][cylinder1._elementsCountAcrossMinor//2+1][1]
        circumferential_element_length = abs(circumferential_element_length)

        armLeftCylinder_basepx = cylinder1._shield.px[upperlimb_elementCountAlong][0][cylinder1._elementsCountAcrossMinor//2]
        armLeftCylinder_basepx[0] += circumferential_element_length

        arm_radius = 0.2
        arm_elements_around = 12

        nx, nd1 = createEllipsePerimeter(armLeftCylinder_basepx, [0,arm_radius,0], [0,0,arm_radius], arm_elements_around, 2*arm_radius)
        nx.pop(-1)
        nd1.pop(-1)

        # square node coordinates.
        arm_elementsCountMinor = 2
        arm_elementsCountMajor = 2
        arm_elementsCountAlong = 1
        px = [[[None for i in range(2*arm_elementsCountMajor+1)] for j in range(2*arm_elementsCountMinor+1)] for k in range(arm_elementsCountAlong+1)]
        pd1 = [[[None for i in range(2*arm_elementsCountMajor+1)] for j in range(2*arm_elementsCountMinor+1)] for k in range(arm_elementsCountAlong+1)]
        pd2 = [[[None for i in range(2*arm_elementsCountMajor+1)] for j in range(2*arm_elementsCountMinor+1)] for k in range(arm_elementsCountAlong+1)]
        pd3 = [[[None for i in range(2*arm_elementsCountMajor+1)] for j in range(2*arm_elementsCountMinor+1)] for k in range(arm_elementsCountAlong+1)]
        nodeId = [[[None for i in range(2*arm_elementsCountMajor+1)] for j in range(2*arm_elementsCountMinor+1)] for k in range(arm_elementsCountAlong+1)]

        d = [circumferential_element_length, 0, 0]
        for n3 in range(arm_elementsCountAlong+1):
            px[n3][0][1] = [nx[2][c] + n3*d[c] for c in range(3)]
            px[n3][0][2] = [nx[3][c] + n3*d[c] for c in range(3)]
            px[n3][0][3] = [nx[4][c] + n3*d[c] for c in range(3)]
            px[n3][1][0] = [nx[1][c] + n3*d[c] for c in range(3)]
            px[n3][1][4] = [nx[5][c] + n3*d[c] for c in range(3)]
            px[n3][2][0] = [nx[0][c] + n3*d[c] for c in range(3)]
            px[n3][2][4] = [nx[6][c] + n3*d[c] for c in range(3)]
            px[n3][3][0] = [nx[11][c] + n3*d[c] for c in range(3)]
            px[n3][3][4] = [nx[7][c] + n3*d[c] for c in range(3)]
            px[n3][4][1] = [nx[10][c] + n3*d[c] for c in range(3)]
            px[n3][4][2] = [nx[9][c] + n3*d[c] for c in range(3)]
            px[n3][4][3] = [nx[8][c] + n3*d[c] for c in range(3)]

            pd1[n3][0][1:4] = nd1[2:5]
            pd1[n3][1][0] = nd1[1]
            pd1[n3][1][4] = nd1[5]
            pd1[n3][2][0] = nd1[0]
            pd1[n3][2][4] = nd1[6]
            pd1[n3][3][0] = nd1[11]
            pd1[n3][3][4] = nd1[7]
            pd1[n3][4][1] = nd1[10]
            pd1[n3][4][2] = nd1[9]
            pd1[n3][4][3] = nd1[8]

            px[n3][1][1] = [armLeftCylinder_basepx[c] + [n3*circumferential_element_length, -2/3*circumferential_element_length, -2/3*circumferential_element_length][c] for c in range(3)]
            px[n3][1][2] = [armLeftCylinder_basepx[c] + [n3*circumferential_element_length, 0, -2/3*circumferential_element_length][c] for c in range(3)]
            px[n3][1][3] = [armLeftCylinder_basepx[c] + [n3*circumferential_element_length, 2/3*circumferential_element_length, -2/3*circumferential_element_length][c] for c in range(3)]
            px[n3][2][1] = [armLeftCylinder_basepx[c] + [n3*circumferential_element_length, -2/3*circumferential_element_length, 0][c] for c in range(3)]
            px[n3][2][2] = [armLeftCylinder_basepx[c] + [n3*circumferential_element_length, 0, 0][c] for c in range(3)]
            px[n3][2][3] = [armLeftCylinder_basepx[c] + [n3*circumferential_element_length, 2/3*circumferential_element_length, 0][c] for c in range(3)]
            px[n3][3][1] = [armLeftCylinder_basepx[c] + [n3*circumferential_element_length, -2/3*circumferential_element_length, 2/3*circumferential_element_length][c] for c in range(3)]
            px[n3][3][2] = [armLeftCylinder_basepx[c] + [n3*circumferential_element_length, 0, 2/3*circumferential_element_length][c] for c in range(3)]
            px[n3][3][3] = [armLeftCylinder_basepx[c] + [n3*circumferential_element_length, 2/3*circumferential_element_length, 2/3*circumferential_element_length][c] for c in range(3)]

            pd1[n3][1][1] = [px[n3][1][2][c] - px[n3][1][1][c] for c in range(3)]
            pd1[n3][1][2] = [px[n3][1][3][c] - px[n3][1][2][c] for c in range(3)]
            pd1[n3][1][3] = [px[n3][1][4][c] - px[n3][1][3][c] for c in range(3)]
            pd1[n3][2][1] = [px[n3][2][2][c] - px[n3][2][1][c] for c in range(3)]
            pd1[n3][2][2] = [px[n3][2][3][c] - px[n3][2][2][c] for c in range(3)]
            pd1[n3][2][3] = [px[n3][2][4][c] - px[n3][2][3][c] for c in range(3)]
            pd1[n3][3][1] = [px[n3][3][2][c] - px[n3][3][1][c] for c in range(3)]
            pd1[n3][3][2] = [px[n3][3][3][c] - px[n3][3][2][c] for c in range(3)]
            pd1[n3][3][3] = [px[n3][3][4][c] - px[n3][3][3][c] for c in range(3)]

            pd3[n3][1][1] = [px[n3][0][1][c] - px[n3][1][1][c] for c in range(3)]
            pd3[n3][1][2] = [px[n3][0][2][c] - px[n3][1][2][c] for c in range(3)]
            pd3[n3][1][3] = [px[n3][0][3][c] - px[n3][1][3][c] for c in range(3)]
            pd3[n3][2][1] = [px[n3][1][1][c] - px[n3][2][1][c] for c in range(3)]
            pd3[n3][2][2] = [px[n3][1][2][c] - px[n3][2][2][c] for c in range(3)]
            pd3[n3][2][3] = [px[n3][1][3][c] - px[n3][2][3][c] for c in range(3)]
            pd3[n3][3][1] = [px[n3][2][1][c] - px[n3][3][1][c] for c in range(3)]
            pd3[n3][3][2] = [px[n3][2][2][c] - px[n3][3][2][c] for c in range(3)]
            pd3[n3][3][3] = [px[n3][2][3][c] - px[n3][3][3][c] for c in range(3)]

        for n3 in range(arm_elementsCountAlong+1):
            for n2 in range(2*arm_elementsCountMinor + 1):
                for n1 in range(2*arm_elementsCountMajor + 1):
                    if px[n3][n2][n1]:
                        pd2[n3][n2][n1] = [circumferential_element_length, 0 ,0]


        for n3 in range(arm_elementsCountAlong+1):
            pd3[n3][0][1] = setMagnitude(crossproduct3(pd1[n3][0][1], pd2[n3][0][1]), circumferential_element_length)
            pd3[n3][0][2] = setMagnitude(crossproduct3(pd1[n3][0][2], pd2[n3][0][2]), circumferential_element_length)
            pd3[n3][0][3] = setMagnitude(crossproduct3(pd1[n3][0][3], pd2[n3][0][3]), circumferential_element_length)
            pd3[n3][1][0] = setMagnitude(crossproduct3(pd1[n3][1][0], pd2[n3][1][0]), circumferential_element_length)
            pd3[n3][1][4] = setMagnitude(crossproduct3(pd1[n3][1][4], pd2[n3][1][4]), circumferential_element_length)
            pd3[n3][2][0] = setMagnitude(crossproduct3(pd1[n3][2][0], pd2[n3][2][0]), circumferential_element_length)
            pd3[n3][2][4] = setMagnitude(crossproduct3(pd1[n3][2][4], pd2[n3][2][4]), circumferential_element_length)
            pd3[n3][3][0] = setMagnitude(crossproduct3(pd1[n3][3][0], pd2[n3][3][0]), circumferential_element_length)
            pd3[n3][3][4] = setMagnitude(crossproduct3(pd1[n3][3][4], pd2[n3][3][4]), circumferential_element_length)
            pd3[n3][4][1] = setMagnitude(crossproduct3(pd1[n3][4][1], pd2[n3][4][1]), circumferential_element_length)
            pd3[n3][4][2] = setMagnitude(crossproduct3(pd1[n3][4][2], pd2[n3][4][2]), circumferential_element_length)
            pd3[n3][4][3] = setMagnitude(crossproduct3(pd1[n3][4][3], pd2[n3][4][3]), circumferential_element_length)

        # generate the nodes
        nodetemplate2 = nodes.createNodetemplate()
        nodetemplate2.defineField(coordinates)
        nodetemplate2.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate2.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate2.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate2.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        cache2 = fieldmodule.createFieldcache()

        for n2 in range(2*arm_elementsCountMinor + 1):
            for n3 in range(arm_elementsCountAlong+1):
                for n1 in range(2*arm_elementsCountMajor + 1):
                    if px[n3][n2][n1]:
                        node = nodes.createNode(nodeIdentifier, nodetemplate2)
                        nodeId[n3][n2][n1] = nodeIdentifier
                        cache2.setNode(node)
                        coordinates.setNodeParameters(cache2, -1, Node.VALUE_LABEL_VALUE, 1, px [n3][n2][n1])
                        coordinates.setNodeParameters(cache2, -1, Node.VALUE_LABEL_D_DS1, 1, pd1[n3][n2][n1])
                        coordinates.setNodeParameters(cache2, -1, Node.VALUE_LABEL_D_DS2, 1, pd2[n3][n2][n1])
                        coordinates.setNodeParameters(cache2, -1, Node.VALUE_LABEL_D_DS3, 1, pd3[n3][n2][n1])
                        nodeIdentifier += 1


        # generate the elements
        elementIdentifier = cylinder1._endElementIdentifier+1
        useCrossDerivatives = False

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftNoCrossDerivatives()
        elementtemplate2 = mesh.createElementtemplate()
        elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate2.defineField(coordinates, -1, eft)

        elementtemplate3 = mesh.createElementtemplate()
        elementtemplate3.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        elementtemplate4 = mesh.createElementtemplate()
        elementtemplate4.setElementShapeType(Element.SHAPE_TYPE_WEDGE12)

        # isEven = (self.elementsCountAcross % 2) == 0
        # e1a = self.elementsCountRim
        # e1b = e1a + 1
        # e1z = self.elementsCountAcross - 1 - self.elementsCountRim
        # e1y = e1z - 1
        # e2a = self.elementsCountRim
        # e2b = self.elementsCountRim + 1
        # e2c = self.elementsCountRim + 2
        # e2z = 2 * self.elementsCountUp - 1 - self.elementsCountRim
        # e2y = e2z - 1
        # e2x = e2z - 2
        wedge= False
        for e3 in range(arm_elementsCountAlong):
            for e2 in range(2*arm_elementsCountMinor):
                for e1 in range(2*arm_elementsCountMajor):
                    eft1 = eft
                    scalefactors = None

                    if (e1>0 and e1<3) and (e2>0 and e2<3):
                        wedge=False
                        nids = [nodeId[e3][e2+1][e1], nodeId[e3][e2 + 1][e1+1], nodeId[e3 + 1][e2+1][e1],
                                nodeId[e3 + 1][e2 + 1][e1+1],
                                nodeId[e3][e2][e1], nodeId[e3][e2][e1 + 1],
                                nodeId[e3 + 1][e2][e1], nodeId[e3 + 1][e2][e1 + 1]]

                    elif (e1 == 0 or e1 == 3) and (e2>0 and e2<3):

                        wedge = False
                        nids = [nodeId[e3][e2 + 1][e1], nodeId[e3][e2 + 1][e1 + 1], nodeId[e3 + 1][e2 + 1][e1],
                                nodeId[e3 + 1][e2 + 1][e1 + 1],
                                nodeId[e3][e2][e1], nodeId[e3][e2][e1 + 1],
                                nodeId[e3 + 1][e2][e1], nodeId[e3 + 1][e2][e1 + 1]]
                        if e1 == 0:
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]
                            remapEftNodeValueLabel(eft1, [1, 3, 5 ,7], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft1, [1, 3, 5 ,7], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS1, [])])
                        if e1 == 3:
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]
                            remapEftNodeValueLabel(eft1, [2, 4, 6 ,8], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [2, 4, 6 ,8], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS3, [])])

                    elif (e2 == 0 or e2 == 3) and (e1>0 and e1<3):

                        wedge = False
                        nids = [nodeId[e3][e2 + 1][e1], nodeId[e3][e2 + 1][e1 + 1], nodeId[e3 + 1][e2 + 1][e1],
                                nodeId[e3 + 1][e2 + 1][e1 + 1],
                                nodeId[e3][e2][e1], nodeId[e3][e2][e1 + 1],
                                nodeId[e3 + 1][e2][e1], nodeId[e3 + 1][e2][e1 + 1]]
                        if e2 == 3:
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]
                            remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS3, [1])])
                    elif e1 == 0 and e2 == 0:
                        wedge = True
                        nids = [nodeId[e3][e2 + 1][e1], nodeId[e3][e2][e1], nodeId[e3 + 1][e2 + 1][e1], nodeId[e3 + 1][e2][e1],
                                nodeId[e3][e2 + 1][e1 + 1], nodeId[e3][e2][e1 + 1], nodeId[e3 + 1][e2 + 1][e1 + 1], nodeId[e3 + 1][e2][e1 + 1]]
                    elif e1 == 3 and e2 == 0:
                        wedge = True
                        nids = [nodeId[e3][e2 + 1][e1 - 1], nodeId[e3][e2 + 1][e1], nodeId[e3 + 1][e2 + 1][e1 - 1], nodeId[e3 + 1][e2 + 1][e1],
                                nodeId[e3][e2 + 1][e1 + 1], nodeId[e3][e2][e1 + 1], nodeId[e3 + 1][e2 + 1][e1 + 1], nodeId[e3 + 1][e2][e1 + 1]]

                    elif e1 == 0 and e2 == 3:
                        wedge = True
                        nids = [nodeId[e3][e2 - 1][e1], nodeId[e3][e2 - 1][e1 + 1], nodeId[e3 + 1][e2 - 1][e1], nodeId[e3 + 1][e2 - 1][e1 + 1],
                                nodeId[e3][e2 - 1][e1 + 1], nodeId[e3][e2][e1 + 1], nodeId[e3 + 1][e2 - 1][e1 + 1], nodeId[e3 + 1][e2][e1 + 1]]
                    elif e1 == 3 and e2 == 3:
                        wedge = True
                        nids = [nodeId[e3][e2 - 1][e1 - 1], nodeId[e3][e2 - 1][e1], nodeId[e3 + 1][e2 - 1][e1 - 1], nodeId[e3 + 1][e2 - 1][e1],
                                nodeId[e3][e2 - 1][e1 - 1], nodeId[e3][e2][e1 - 1], nodeId[e3 + 1][e2 - 1][e1 - 1], nodeId[e3 + 1][e2][e1 - 1]]



                    # if (e2 < e2b) or (e2 > e2y):
                    #     if (e1 < e1b) or (e1 > e1y):
                    #         continue  # no element due to triple point closure
                    #     if (e2 < e2a) or (e2 > e2z):
                    #         if self._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                    #             if e2 < e2a:
                    #                 nids = [self.nodeId[e3][e2 + 1][e1], self.nodeId[e3][e2 + 1][e1 + 1],
                    #                         self.nodeId[e3 + 1][e2 + 1][e1], self.nodeId[e3 + 1][e2 + 1][e1 + 1],
                    #                         self.nodeId[e3][e2][e1], self.nodeId[e3][e2][e1 + 1],
                    #                         self.nodeId[e3 + 1][e2][e1], self.nodeId[e3 + 1][e2][e1 + 1]]
                    #             elif e2 > e2z:
                    #                 nids = [self.nodeId[e3][e2][e1 + 1], self.nodeId[e3][e2][e1],
                    #                         self.nodeId[e3 + 1][e2][e1 + 1], self.nodeId[e3 + 1][e2][e1],
                    #                         self.nodeId[e3][e2 + 1][e1 + 1], self.nodeId[e3][e2 + 1][e1],
                    #                         self.nodeId[e3 + 1][e2 + 1][e1 + 1], self.nodeId[e3 + 1][e2 + 1][e1]]
                    #     elif (e2 == e2a) or (e2 == e2z):
                    #         # bottom and top row elements
                    #         if self._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                    #             if e2 == e2a:
                    #                 eft1 = tricubichermite.createEftNoCrossDerivatives()
                    #                 setEftScaleFactorIds(eft1, [1], [])
                    #                 scalefactors = [-1.0]
                    #                 remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                    #                                        [(Node.VALUE_LABEL_D_DS3, [1])])
                    #                 remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS3,
                    #                                        [(Node.VALUE_LABEL_D_DS1, [])])

                    if wedge:
                        # elementtemplate4.defineField(coordinates, -1, eft)
                        # element = mesh.createElement(elementIdentifier, elementtemplate4)
                        pass
                    else:
                        if eft1 is not eft:
                            elementtemplate3.defineField(coordinates, -1, eft1)
                            element = mesh.createElement(elementIdentifier, elementtemplate3)
                        else:
                            element = mesh.createElement(elementIdentifier, elementtemplate2)
                        result2 = element.setNodesByIdentifier(eft1, nids)
                    if scalefactors:
                        result3 = element.setScaleFactors(eft1, scalefactors)
                    else:
                        result3 = 7

                    elementIdentifier += 1


        # for meshGroup in meshGroups:
        #     meshGroup.addElement(element)
        for n2 in range(2*arm_elementsCountMinor + 1):
            for n1 in range(2*arm_elementsCountMajor + 1):
                nodeId[1][n2][n1] = nodeId[0][n2][n1]

        nodeId[0][0][1:3] = [80,81,82]
        nodeId[0][1][0:5] = [86, 87, 88, 89, 90]
        nodeId[0][2][0:5] = [93, 94, 95, 96, 97]
        nodeId[0][3][0:5] = [100, 101, 102, 103, 104]
        nodeId[0][4][1:3] = [108, 109, 110]

        wedge = False
        for e3 in range(arm_elementsCountAlong):
            for e2 in range(2 * arm_elementsCountMinor):
                for e1 in range(2 * arm_elementsCountMajor):
                    eft1 = eft
                    scalefactors = None

                    if (e1 > 0 and e1 < 3) and (e2 > 0 and e2 < 3):
                        wedge = False
                        nids = [nodeId[e3][e2 + 1][e1], nodeId[e3][e2 + 1][e1 + 1],
                                nodeId[e3 + 1][e2 + 1][e1],
                                nodeId[e3 + 1][e2 + 1][e1 + 1],
                                nodeId[e3][e2][e1], nodeId[e3][e2][e1 + 1],
                                nodeId[e3 + 1][e2][e1], nodeId[e3 + 1][e2][e1 + 1]]
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [-1.0]
                        remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2,
                                               [(Node.VALUE_LABEL_D_DS3, [])])
                        remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS3,
                                               [(Node.VALUE_LABEL_D_DS2, [1])])

                    elif (e1 == 0 or e1 == 3) and (e2 > 0 and e2 < 3):

                        wedge = False
                        nids = [nodeId[e3][e2 + 1][e1], nodeId[e3][e2 + 1][e1 + 1],
                                nodeId[e3 + 1][e2 + 1][e1],
                                nodeId[e3 + 1][e2 + 1][e1 + 1],
                                nodeId[e3][e2][e1], nodeId[e3][e2][e1 + 1],
                                nodeId[e3 + 1][e2][e1], nodeId[e3 + 1][e2][e1 + 1]]
                        if e1 == 0:
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS1, [])])

                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [])])
                        if e1 == 3:
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS3, [])])

                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])


                    elif (e2 == 0 or e2 == 3) and (e1 > 0 and e1 < 3):

                        wedge = False
                        nids = [nodeId[e3][e2 + 1][e1], nodeId[e3][e2 + 1][e1 + 1],
                                nodeId[e3 + 1][e2 + 1][e1],
                                nodeId[e3 + 1][e2 + 1][e1 + 1],
                                nodeId[e3][e2][e1], nodeId[e3][e2][e1 + 1],
                                nodeId[e3 + 1][e2][e1], nodeId[e3 + 1][e2][e1 + 1]]

                        if e2 == 0:
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]

                            remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS3, [])])

                        if e2 == 3:
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]
                            remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS3, [1])])

                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])

                    elif e1 == 0 and e2 == 0:
                        wedge = True
                        nids = [nodeId[e3][e2 + 1][e1], nodeId[e3][e2][e1], nodeId[e3 + 1][e2 + 1][e1],
                                nodeId[e3 + 1][e2][e1],
                                nodeId[e3][e2 + 1][e1 + 1], nodeId[e3][e2][e1 + 1],
                                nodeId[e3 + 1][e2 + 1][e1 + 1], nodeId[e3 + 1][e2][e1 + 1]]
                    elif e1 == 3 and e2 == 0:
                        wedge = True
                        nids = [nodeId[e3][e2 + 1][e1 - 1], nodeId[e3][e2 + 1][e1],
                                nodeId[e3 + 1][e2 + 1][e1 - 1], nodeId[e3 + 1][e2 + 1][e1],
                                nodeId[e3][e2 + 1][e1 + 1], nodeId[e3][e2][e1 + 1],
                                nodeId[e3 + 1][e2 + 1][e1 + 1], nodeId[e3 + 1][e2][e1 + 1]]

                    elif e1 == 0 and e2 == 3:
                        wedge = True
                        nids = [nodeId[e3][e2 - 1][e1], nodeId[e3][e2 - 1][e1 + 1],
                                nodeId[e3 + 1][e2 - 1][e1], nodeId[e3 + 1][e2 - 1][e1 + 1],
                                nodeId[e3][e2 - 1][e1 + 1], nodeId[e3][e2][e1 + 1],
                                nodeId[e3 + 1][e2 - 1][e1 + 1], nodeId[e3 + 1][e2][e1 + 1]]
                    elif e1 == 3 and e2 == 3:
                        wedge = True
                        nids = [nodeId[e3][e2 - 1][e1 - 1], nodeId[e3][e2 - 1][e1],
                                nodeId[e3 + 1][e2 - 1][e1 - 1], nodeId[e3 + 1][e2 - 1][e1],
                                nodeId[e3][e2 - 1][e1 - 1], nodeId[e3][e2][e1 - 1],
                                nodeId[e3 + 1][e2 - 1][e1 - 1], nodeId[e3 + 1][e2][e1 - 1]]
                    if wedge:
                        # elementtemplate4.defineField(coordinates, -1, eft)
                        # element = mesh.createElement(elementIdentifier, elementtemplate4)
                        pass
                    else:
                        if eft1 is not eft:
                            elementtemplate3.defineField(coordinates, -1, eft1)
                            element = mesh.createElement(elementIdentifier, elementtemplate3)
                        else:
                            element = mesh.createElement(elementIdentifier, elementtemplate2)
                        result2 = element.setNodesByIdentifier(eft1, nids)
                    if scalefactors:
                        result3 = element.setScaleFactors(eft1, scalefactors)
                    else:
                        result3 = 7

                    elementIdentifier += 1







        a = 1

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCountAcrossMajor = options['Refine number of elements across major']
        refineElementsCountAlong = options['Refine number of elements along']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCountAcrossMajor, refineElementsCountAlong, refineElementsCountAcrossMajor)

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
        mesh2d = fm.findMeshByDimension(2)

        bodyGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("body"))
        coreGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("core"))
        non_coreGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("non core"))
        abdomenGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("abdomen"))
        thoraxGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("thorax"))
        neckGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("neck core"))

        skinGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("skin epidermis"))
        coreBoundaryGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("core boundary"))
        diaphragmGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("diaphragm"))
        spinalCordGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("spinal cord"))

        is_exterior = fm.createFieldIsExterior()
        is_on_face_xi3_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1)
        is_skin = fm.createFieldAnd(is_exterior, is_on_face_xi3_1)

        skinMeshGroup = skinGroup.getMeshGroup(mesh2d)
        skinMeshGroup.addElementsConditional(is_skin)

        is_core_boundary = fm.createFieldAnd(coreGroup.getGroup(), non_coreGroup.getGroup())
        coreBoundaryMeshGroup = coreBoundaryGroup.getMeshGroup(mesh2d)
        coreBoundaryMeshGroup.addElementsConditional(is_core_boundary)

        is_diaphragm = fm.createFieldAnd(abdomenGroup.getGroup(), thoraxGroup.getGroup())
        diaphragmMeshGroup = diaphragmGroup.getMeshGroup(mesh2d)
        diaphragmMeshGroup.addElementsConditional(is_diaphragm)

        # spinal cord
        coordinates = fm.findFieldByName('coordinates').castFiniteElement()
        zero = fm.createFieldConstant(0)
        zero_m = fm.createFieldConstant(-0.01)
        zero_p = fm.createFieldConstant(0.01)
        comp2 = cls.axis2.index(max(cls.axis2)) + 1
        ax2_comp = fm.createFieldComponent(coordinates, comp2)
        ax2_gt_zero_m = fm.createFieldGreaterThan(ax2_comp, zero_m)
        ax2_lt_zero_p = fm.createFieldLessThan(ax2_comp, zero_p)
        ax2_gt_zero_xi10 = fm.createFieldAnd(fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0), ax2_gt_zero_m)
        ax2_lt_zero_xi10 = fm.createFieldAnd(fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0), ax2_lt_zero_p)
        is_ax2_zero = fm.createFieldAnd(ax2_lt_zero_xi10, ax2_gt_zero_xi10)
        comp3 = cls.axis3.index(max(cls.axis3)) + 1
        ax3_comp = fm.createFieldComponent(coordinates, comp3)
        ax3_positive = fm.createFieldGreaterThan(ax3_comp, zero)
        is_ax2_zero_ax3_positive = fm.createFieldAnd(is_ax2_zero, ax3_positive)
        is_abdomen_thorax = fm.createFieldAdd(abdomenGroup.getGroup(), thoraxGroup.getGroup())
        is_abdomen_thorax_neck = fm.createFieldAdd(is_abdomen_thorax, neckGroup.getGroup())
        is_abdomen_thorax_neck_boundary = fm.createFieldAnd(is_core_boundary, is_abdomen_thorax_neck)
        is_spinal_cord = fm.createFieldAnd(is_ax2_zero_ax3_positive, is_abdomen_thorax_neck_boundary)

        mesh1d = fm.findMeshByDimension(1)
        spinalCordMeshGroup = spinalCordGroup.getMeshGroup(mesh1d)
        spinalCordMeshGroup.addElementsConditional(is_spinal_cord)
