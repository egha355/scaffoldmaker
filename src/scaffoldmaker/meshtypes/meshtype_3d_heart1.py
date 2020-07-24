"""
Generates a 3-D heart model including ventricles, base and atria.
"""

from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldNodeGroup, findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.utils.zinc.finiteelement import getMaximumElementIdentifier, getMaximumNodeIdentifier
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm, mergeAnnotationGroups
from scaffoldmaker.annotation.heart_terms import get_heart_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.meshtypes.meshtype_3d_heartatria1 import MeshType_3d_heartatria1
from scaffoldmaker.meshtypes.meshtype_3d_heartventriclesbase1 import MeshType_3d_heartventriclesbase1
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.meshrefinement import MeshRefinement

class MeshType_3d_heart1(Scaffold_base):
    '''
    Generates a 3-D heart model including ventricles, base and atria.
    '''

    @staticmethod
    def getName():
        return '3D Heart 1'

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
        options = MeshType_3d_heartatria1.getDefaultOptions(parameterSetName)
        # ventricles overrides some atria options, so update from it:
        optionsVentricles = MeshType_3d_heartventriclesbase1.getDefaultOptions(parameterSetName)
        options.update(optionsVentricles)
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = MeshType_3d_heartventriclesbase1.getOrderedOptionNames()
        optionNamesAtria = MeshType_3d_heartatria1.getOrderedOptionNames()
        # insert new numbers of elements from atria; others came with ventriclesbase1
        optionNames.insert(7, 'Number of elements along vena cava inlet')
        optionNames.insert(8, 'Number of elements over atria')
        # remove number of elements, unit scale and dependent options from atria options
        for key in [
            'Number of elements along vena cava inlet',
            'Number of elements around atrial septum',
            'Number of elements around left atrium free wall',
            'Number of elements around right atrium free wall',
            'Number of elements over atria',
            'Unit scale',
            'Aorta outer plus diameter']:
            optionNamesAtria.remove(key)
        # remove duplicates from main options
        for optionName in optionNamesAtria:
            if optionName in optionNames:
                optionNames.remove(optionName)
        # add atria options all together
        optionNames += optionNamesAtria
        # want refinement options last
        for optionName in [
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through LV wall',
            'Refine number of elements through wall']:
            optionNames.remove(optionName)
            optionNames.append(optionName)
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        return MeshType_3d_heartatria1.getOptionValidScaffoldTypes(optionName)

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        return MeshType_3d_heartatria1.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType)

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        return MeshType_3d_heartatria1.getOptionScaffoldPackage(optionName, scaffoldType, parameterSetName)

    @staticmethod
    def checkOptions(options):
        '''
        :return:  True if dependent options changed, otherwise False. This
        happens where two or more options must change together to be valid.
        '''
        dependentChanges = MeshType_3d_heartventriclesbase1.checkOptions(options) \
            or MeshType_3d_heartatria1.checkOptions(options)
        # set dependent outer diameter used in atria2
        options['Aorta outer plus diameter'] = options['LV outlet inner diameter'] + 2.0*options['LV outlet wall thickness']
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        # set dependent outer diameter used in atria2
        options['Aorta outer plus diameter'] = options['LV outlet inner diameter'] + 2.0*options['LV outlet wall thickness']
        elementsCountAroundAtrialSeptum = options['Number of elements around atrial septum']
        elementsCountAroundLeftAtriumFreeWall = options['Number of elements around left atrium free wall']
        elementsCountAroundLeftAtrium = elementsCountAroundLeftAtriumFreeWall + elementsCountAroundAtrialSeptum
        elementsCountAroundRightAtriumFreeWall = options['Number of elements around right atrium free wall']
        elementsCountAroundRightAtrium = elementsCountAroundRightAtriumFreeWall + elementsCountAroundAtrialSeptum
        useCrossDerivatives = False

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)
        cache = fm.createFieldcache()

        mesh = fm.findMeshByDimension(3)

        # generate heartventriclesbase1 model and put atria1 on it
        ventriclesAnnotationGroups = MeshType_3d_heartventriclesbase1.generateBaseMesh(region, options)
        atriaAnnotationGroups = MeshType_3d_heartatria1.generateBaseMesh(region, options)
        annotationGroups = mergeAnnotationGroups(ventriclesAnnotationGroups, atriaAnnotationGroups)
        lFibrousRingGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_heart_term("left fibrous ring"))
        rFibrousRingGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_heart_term("right fibrous ring"))

        # annotation fiducial points
        markerGroup = findOrCreateFieldGroup(fm, "marker")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)

        ##############
        # Create nodes
        ##############

        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)

        # discover left and right fibrous ring nodes from ventricles and atria
        # because nodes are iterated in identifier order, the lowest and first are on the lv outlet cfb, right and left on lower outer layers
        # left fibrous ring
        lavNodeId = [ [ [], [] ], [ [], [] ] ]  # [n3][n2][n1]
        iter = lFibrousRingGroup.getNodesetGroup(nodes).createNodeiterator()
        # left fibrous ring, bottom row
        cfbNodeId = iter.next().getIdentifier()
        cfbLeftNodeId = iter.next().getIdentifier()
        for n1 in range(elementsCountAroundLeftAtrium):
            lavNodeId[0][0].append(iter.next().getIdentifier())
        lavNodeId[1][0].append(cfbNodeId)
        lavNodeId[1][0].append(cfbLeftNodeId)
        for n1 in range(elementsCountAroundLeftAtriumFreeWall - 1):
            lavNodeId[1][0].append(iter.next().getIdentifier())
        for n1 in range(elementsCountAroundAtrialSeptum - 1):
            lavNodeId[1][0].append(None)  # no outer node on interatrial septum
        # left fibrous ring, top row
        for n1 in range(elementsCountAroundLeftAtrium):
            lavNodeId[0][1].append(iter.next().getIdentifier())
        for n1 in range(elementsCountAroundLeftAtriumFreeWall + 1):
            lavNodeId[1][1].append(iter.next().getIdentifier())
        for n1 in range(elementsCountAroundAtrialSeptum - 1):
            lavNodeId[1][1].append(None)  # no outer node on interatrial septum
        # right fibrous ring
        ravNodeId = [ [ [], [] ], [ [], [] ] ]  # [n3][n2][n1]
        iter = rFibrousRingGroup.getNodesetGroup(nodes).createNodeiterator()
        cfbNodeId = iter.next().getIdentifier()
        cfbRightNodeId = iter.next().getIdentifier()
        # right fibrous ring, bottom row
        for n1 in range(elementsCountAroundRightAtrium):
            ravNodeId[0][0].append(iter.next().getIdentifier())
        for n1 in range(elementsCountAroundRightAtriumFreeWall - 1):
            ravNodeId[1][0].append(iter.next().getIdentifier())
        ravNodeId[1][0].append(cfbRightNodeId)
        ravNodeId[1][0].append(cfbNodeId)
        for n1 in range(elementsCountAroundAtrialSeptum - 1):
            ravNodeId[1][0].append(None)  # no outer node on interatrial septum
        # right fibrous ring, top row
        for n1 in range(elementsCountAroundRightAtrium):
            ravNodeId[0][1].append(iter.next().getIdentifier())
        cfbUpperNodeId = iter.next().getIdentifier()  # cfb from left will be first
        for n1 in range(elementsCountAroundRightAtriumFreeWall):
            ravNodeId[1][1].append(iter.next().getIdentifier())
        ravNodeId[1][1].append(cfbUpperNodeId)
        for n1 in range(elementsCountAroundAtrialSeptum - 1):
            ravNodeId[1][1].append(None)  # no outer node on interatrial septum

        #for n2 in range(2):
        #    print('n2', n2)
        #    print('lavNodeId[0]', lavNodeId[0][n2])
        #    print('lavNodeId[1]', lavNodeId[1][n2])
        #    print('ravNodeId[0]', ravNodeId[0][n2])
        #    print('ravNodeId[1]', ravNodeId[1][n2])

        #################
        # Create elements
        #################

        lFibrousRingMeshGroup = lFibrousRingGroup.getMeshGroup(mesh)
        rFibrousRingMeshGroup = rFibrousRingGroup.getMeshGroup(mesh)

        elementIdentifier = getMaximumElementIdentifier(mesh) + 1

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # create fibrous ring elements

        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives, linearAxis = 2, d_ds1 = Node.VALUE_LABEL_D_DS1, d_ds2 = Node.VALUE_LABEL_D_DS3)
        eftFibrousRing = bicubichermitelinear.createEftBasic()

        # left fibrous ring, starting at crux / collapsed posterior interatrial sulcus
        cruxElementId = None
        for e in range(-1, elementsCountAroundLeftAtriumFreeWall):
            eft1 = eftFibrousRing
            n1 = e
            nids = [
                lavNodeId[0][0][n1], lavNodeId[0][0][n1 + 1], lavNodeId[0][1][n1], lavNodeId[0][1][n1 + 1],
                lavNodeId[1][0][n1], lavNodeId[1][0][n1 + 1], lavNodeId[1][1][n1], lavNodeId[1][1][n1 + 1]]
            scalefactors = None
            meshGroups = [ lFibrousRingMeshGroup ]

            if e == -1:
                # interatrial groove straddles left and right atria, collapsed to 6 node wedge
                nids[0] = ravNodeId[0][0][elementsCountAroundRightAtriumFreeWall]
                nids[2] = ravNodeId[0][1][elementsCountAroundRightAtriumFreeWall]
                nids.pop(6)
                nids.pop(4)
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                # reverse d3 on cfb:
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [0] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                ln_map = [ 1, 2, 3, 4, 5, 5, 6, 6 ]
                remapEftLocalNodes(eft1, 6, ln_map)
                meshGroups += [ rFibrousRingMeshGroup ]
            elif e == 0:
                # general linear map d3 adjacent to collapsed sulcus
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                # reverse d1, d3 on cfb, left cfb:
                scaleEftNodeValueLabels(eft1, [ 6 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
            elif e == 1:
                # reverse d1, d3 on left cfb:
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1]) ])
            elif e == (elementsCountAroundLeftAtriumFreeWall - 1):
                # general linear map d3 adjacent to collapsed sulcus
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])

            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            result3 = element.setScaleFactors(eft1, scalefactors) if scalefactors else None
            #print('create element fibrous ring left', elementIdentifier, result, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # right fibrous ring, starting at crux / collapsed posterior interatrial sulcus
        for e in range(-1, elementsCountAroundRightAtriumFreeWall):
            eft1 = eftFibrousRing
            n1 = e
            nids = [
                ravNodeId[0][0][n1], ravNodeId[0][0][n1 + 1], ravNodeId[0][1][n1], ravNodeId[0][1][n1 + 1],
                ravNodeId[1][0][n1], ravNodeId[1][0][n1 + 1], ravNodeId[1][1][n1], ravNodeId[1][1][n1 + 1]]
            scalefactors = None
            meshGroups = [ rFibrousRingMeshGroup ]

            if e == -1:
                # interatrial groove straddles left and right atria, collapsed to 6 node wedge
                nids[0] = lavNodeId[0][0][elementsCountAroundLeftAtriumFreeWall]
                nids[2] = lavNodeId[0][1][elementsCountAroundLeftAtriumFreeWall]
                nids.pop(6)
                nids.pop(4)
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                ln_map = [ 1, 2, 3, 4, 5, 5, 6, 6 ]
                remapEftLocalNodes(eft1, 6, ln_map)
                meshGroups += [ lFibrousRingMeshGroup ]
                cruxElementId = elementIdentifier
            elif e == 0:
                # general linear map d3 adjacent to collapsed crux/posterior sulcus
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
            elif e == (elementsCountAroundRightAtriumFreeWall - 2):
                # reverse d1, d3 on right cfb:
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1]) ])
            elif e == (elementsCountAroundRightAtriumFreeWall - 1):
                # general linear map d3 adjacent to collapsed cfb/anterior sulcus
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                # reverse d1, d3 on right cfb, cfb:
                scaleEftNodeValueLabels(eft1, [ 5 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])

            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            result3 = element.setScaleFactors(eft1, scalefactors) if scalefactors else None
            #print('create element fibrous ring right', elementIdentifier, result, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # fibrous ring septum:
        meshGroups = [ lFibrousRingMeshGroup, rFibrousRingMeshGroup ]
        for e in range(elementsCountAroundAtrialSeptum):
            eft1 = eftFibrousRing
            nlm = e - elementsCountAroundAtrialSeptum
            nlp = nlm + 1
            nrm = -e
            nrp = nrm - 1
            nids = [ lavNodeId[0][0][nlm], lavNodeId[0][0][nlp], lavNodeId[0][1][nlm], lavNodeId[0][1][nlp],
                     ravNodeId[0][0][nrm], ravNodeId[0][0][nrp], ravNodeId[0][1][nrm], ravNodeId[0][1][nrp] ]

            eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
            setEftScaleFactorIds(eft1, [1], [])
            scalefactors = [ -1.0 ]
            if e == 0:
                # general linear map d3 adjacent to collapsed posterior interventricular sulcus
                scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                scaleEftNodeValueLabels(eft1, [ 6, 8 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
            elif e == (elementsCountAroundAtrialSeptum - 1):
                # general linear map d3 adjacent to cfb
                scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                scaleEftNodeValueLabels(eft1, [ 5, 7 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
            else:
                scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])

            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            result3 = element.setScaleFactors(eft1, scalefactors) if scalefactors else None
            #print('create element fibrous ring septum', elementIdentifier, result, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # annotation fiducial points
        cruxElement = mesh.findElementByIdentifier(cruxElementId)
        cruxXi = [ 0.5, 0.5, 1.0 ]
        cache.setMeshLocation(cruxElement, cruxXi)
        result, cruxCoordinates = coordinates.evaluateReal(cache, 3)
        markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
        nodeIdentifier += 1
        cache.setNode(markerPoint)
        markerName.assignString(cache, "crux of heart")
        markerLocation.assignMeshLocation(cache, cruxElement, cruxXi)

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        elementsCountAroundAtrialSeptum = options['Number of elements around atrial septum']
        elementsCountAroundLeftAtriumFreeWall = options['Number of elements around left atrium free wall']
        elementsCountAroundRightAtriumFreeWall = options['Number of elements around right atrium free wall']
        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughWall = options['Refine number of elements through wall']
        MeshType_3d_heartventriclesbase1.refineMesh(meshrefinement, options)
        MeshType_3d_heartatria1.refineMesh(meshrefinement, options)
        element = meshrefinement._sourceElementiterator.next()
        startFibrousRingElementIdentifier = element.getIdentifier()
        lastFibrousRingElementIdentifier = startFibrousRingElementIdentifier + elementsCountAroundLeftAtriumFreeWall + elementsCountAroundRightAtriumFreeWall + elementsCountAroundAtrialSeptum + 1
        i = -1
        while element.isValid():
            numberInXi1 = refineElementsCountSurface
            numberInXi2 = 1  # since fibrous ring is thin
            numberInXi3 = refineElementsCountThroughWall
            elementIdentifier = element.getIdentifier()
            if i in [ -1, elementsCountAroundLeftAtriumFreeWall ]:
                numberInXi1 = refineElementsCountThroughWall
            meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
            if elementIdentifier == lastFibrousRingElementIdentifier:
                return  # finish on last so can continue elsewhere
            element = meshrefinement._sourceElementiterator.next()
            i += 1

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
        # add epicardium of fibrous ring
        fm = region.getFieldmodule()
        MeshType_3d_heartventriclesbase1.defineFaceAnnotations(region, options, annotationGroups)
        MeshType_3d_heartatria1.defineFaceAnnotations(region, options, annotationGroups)
        lFibrousRingGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("left fibrous ring"))
        rFibrousRingGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("right fibrous ring"))
        epiGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("epicardium"))
        mesh2d = fm.findMeshByDimension(2)
        is_exterior_face_xi3_1 = fm.createFieldAnd(fm.createFieldIsExterior(), fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_non_septal_fibrous_ring = fm.createFieldXor(lFibrousRingGroup.getFieldElementGroup(mesh2d), rFibrousRingGroup.getFieldElementGroup(mesh2d))
        is_non_septal_fibrous_ring_epi = fm.createFieldAnd(is_non_septal_fibrous_ring, is_exterior_face_xi3_1)
        epiGroup.getMeshGroup(mesh2d).addElementsConditional(is_non_septal_fibrous_ring_epi)
