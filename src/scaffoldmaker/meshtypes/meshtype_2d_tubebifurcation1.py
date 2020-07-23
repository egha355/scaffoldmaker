"""
Generates a 2-D tube bifurcation mesh.
"""

from __future__ import division
import math
from opencmiss.utils.maths.vectorops import cross, mult, normalize, sub
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.bifurcation import get_tube_bifurcation_connection_elements_counts, \
    make_tube_bifurcation_points, make_tube_bifurcation_elements_2d
from scaffoldmaker.utils.geometry import createCirclePoints


class MeshType_2d_tubebifurcation1(Scaffold_base):
    '''
    Generates a 2-D tube bifurcation mesh.
    '''

    @staticmethod
    def getName():
        return '2D Tube Bifurcation 1'

    @staticmethod
    def getParameterSetNames():
        return ['Default']

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        options = {}
        options['Number of elements around parent'] = 6
        options['Number of elements around child 1'] = 6
        options['Number of elements around child 2'] = 6
        options['Parent diameter'] = 0.4
        options['Parent length'] = 0.5
        options['Parent length scale factor'] = 0.8
        options['Child 1 angle degrees'] = 45.0
        options['Child 1 diameter'] = 0.36
        options['Child 1 length'] = 0.5
        options['Child 1 length scale factor'] = 0.8
        options['Child 2 angle degrees'] = -45.0
        options['Child 2 diameter'] = 0.2
        options['Child 2 length'] = 0.5
        options['Child 2 length scale factor'] = 0.8
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around parent',
            'Number of elements around child 1',
            'Number of elements around child 2',
            'Parent diameter',
            'Parent length',
            'Parent length scale factor',
            'Child 1 angle degrees',
            'Child 1 diameter',
            'Child 1 length',
            'Child 1 length scale factor',
            'Child 2 angle degrees', 
            'Child 2 diameter',
            'Child 2 length',
            'Child 2 length scale factor']

    @staticmethod
    def checkOptions(options):
        '''
        :return:  True if dependent options changed, otherwise False. This
        happens where two or more options must change together to be valid.
        '''
        dependentChanges = False
        for key in [
            'Number of elements around parent',
            'Number of elements around child 1',
            'Number of elements around child 2']:
            if options[key] < 4:
                options[key] = 4
        paCount = options['Number of elements around parent']
        c1Count = options['Number of elements around child 1']
        c2Count = options['Number of elements around child 2']
        if (paCount + c1Count + c2Count) % 2:
            c2Count += 1
        pac1Count, pac2Count, c1c2Count = get_tube_bifurcation_connection_elements_counts(paCount, c1Count, c2Count)
        if pac1Count < 2:
            c2Count -= 2*(2 - pac1Count)
        elif pac2Count < 2:
            c2Count += 2*(2 - pac2Count)
        elif c1c2Count < 2:
            c2Count += 2*(2 - c1c2Count)
        if c2Count != options['Number of elements around child 2']:
            options['Number of elements around child 2'] = c2Count
            dependentChanges = True
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base bicubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        paCount = options['Number of elements around parent']
        c1Count = options['Number of elements around child 1']
        c2Count = options['Number of elements around child 2']
        parentRadius = 0.5*options['Parent diameter']
        parentLength = options['Parent length']
        parentLengthScaleFactor = options['Parent length scale factor']
        child1AngleRadians = math.radians(options['Child 1 angle degrees'])
        child1Radius = 0.5*options['Child 1 diameter']
        child1Length = options['Child 1 length']
        child1LengthScaleFactor = options['Child 1 length scale factor']
        child2AngleRadians = math.radians(options['Child 2 angle degrees'])
        child2Radius = 0.5*options['Child 2 diameter']
        child2Length = options['Child 2 length']
        child2LengthScaleFactor = options['Child 2 length scale factor']
        useCrossDerivatives = False

        # centres:
        paCentre = [ 0.0, 0.0, -parentLength ]
        c1Centre = [ child1Length*math.sin(child1AngleRadians), 0.0, child1Length*math.cos(child1AngleRadians) ]
        c2Centre = [ child2Length*math.sin(child2AngleRadians), 0.0, child2Length*math.cos(child2AngleRadians) ]
        c12 = sub(c1Centre, c2Centre)

        pac1Count, pac2Count, c1c2Count = get_tube_bifurcation_connection_elements_counts(paCount, c1Count, c2Count)

        # parent ring
        paAxis3 = [ 0.0, 0.0, parentLength ]
        paAxis2 = mult(normalize(cross(paAxis3, c12)), parentRadius)
        paAxis1 = mult(normalize(cross(paAxis2, paAxis3)), parentRadius)
        paStartRadians = -math.pi*(pac1Count/paCount)
        pax, pad1 = createCirclePoints(paCentre, paAxis1, paAxis2, paCount, paStartRadians)
        pad2 = [ mult(paAxis3, parentLengthScaleFactor) ]*paCount
        # child 1 ring
        c1Axis3 = c1Centre
        c1Axis2 = mult(normalize(cross(c1Axis3, c12)), child1Radius)
        c1Axis1 = mult(normalize(cross(c1Axis2, c1Axis3)), child1Radius)
        c1StartRadians = -math.pi*(pac1Count/c1Count)
        c1x, c1d1 = createCirclePoints(c1Centre, c1Axis1, c1Axis2, c1Count, c1StartRadians)
        c1d2 = [  mult(c1Axis3, child1LengthScaleFactor) ]*c1Count
        # child 2 ring
        c2Axis3 = c2Centre
        c2Axis2 = mult(normalize(cross(c2Axis3, c12)), child2Radius)
        c2Axis1 = mult(normalize(cross(c2Axis2, c2Axis3)), child2Radius)
        c2StartRadians = -math.pi*(c1c2Count/c2Count)
        c2x, c2d1 = createCirclePoints(c2Centre, c2Axis1, c2Axis2, c2Count, c2StartRadians)
        c2d2 = [ mult(c2Axis3, child2LengthScaleFactor) ]*c2Count

        rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex = \
            make_tube_bifurcation_points(paCentre, pax, pad2, c1Centre, c1x, c1d2, c2Centre, c2x, c2d2)

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)
        cache = fm.createFieldcache()

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        ##############
        # Create nodes
        ##############

        nodeIdentifier = 1

        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        paNodeId = []
        for n in range(paCount):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, pax [n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, pad1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, pad2[n])
            paNodeId.append(nodeIdentifier)
            nodeIdentifier = nodeIdentifier + 1
        roNodeId = []
        for n in range(len(rox)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rox [n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rod1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rod2[n])
            roNodeId.append(nodeIdentifier)
            nodeIdentifier = nodeIdentifier + 1
        coNodeId = []
        for n in range(len(cox)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cox [n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cod1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cod2[n])
            coNodeId.append(nodeIdentifier)
            nodeIdentifier = nodeIdentifier + 1
        c1NodeId = []
        for n in range(c1Count):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, c1x [n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, c1d1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, c1d2[n])
            c1NodeId.append(nodeIdentifier)
            nodeIdentifier = nodeIdentifier + 1
        c2NodeId = []
        for n in range(c2Count):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, c2x [n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, c2d1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, c2d2[n])
            c2NodeId.append(nodeIdentifier)
            nodeIdentifier = nodeIdentifier + 1

        #################
        # Create elements
        #################

        elementIdentifier = 1
        elementIdentifier = make_tube_bifurcation_elements_2d(region, coordinates, elementIdentifier,
            paNodeId, paStartIndex, c1NodeId, c1StartIndex, c2NodeId, c2StartIndex, roNodeId, coNodeId, useCrossDerivatives)

        return []
