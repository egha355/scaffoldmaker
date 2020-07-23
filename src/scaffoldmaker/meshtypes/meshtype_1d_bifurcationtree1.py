"""
Generates a 1-D tree of bifurcating curves with radius.
"""

from __future__ import division
from math import cos, radians, sin, sqrt
from opencmiss.utils.maths.vectorops import add, cross, magnitude, mult, normalize, sub
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldFiniteElement
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base


class MeshType_1d_bifurcationtree1(Scaffold_base):
    '''
    Generates a 1-D tree of bifurcating curves with radius.
    '''

    @staticmethod
    def getName():
        return '1D Bifurcation Tree 1'

    @staticmethod
    def getParameterSetNames():
        return ['Default']

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        options = {}
        options['Number of generations'] = 8
        options['Root diameter'] = 0.25
        options['Root length'] = 1.0
        options['Fork angle degrees'] = 10.0
        options['Fork diameter ratio'] = 0.9
        options['Branch arc angle degrees'] = 40.0
        options['Branch diameter ratio'] = 0.8
        options['Branch length ratio'] = 0.8
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of generations',
            'Root diameter',
            'Root length',
            'Fork angle degrees',
            'Fork diameter ratio',
            'Branch arc angle degrees',
            'Branch diameter ratio',
            'Branch length ratio']

    @staticmethod
    def checkOptions(options):
        '''
        :return:  True if dependent options changed, otherwise False. This
        happens where two or more options must change together to be valid.
        '''
        dependentChanges = False
        for key in ['Number of generations']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Root diameter',
            'Root length',
            'Fork angle degrees',
            'Fork diameter ratio',
            'Branch arc angle degrees',
            'Branch diameter ratio',
            'Branch length ratio']:
            if options[key] < 0.0:
                options[key] = 0.0
        return dependentChanges

    @classmethod
    def generateBifurcationTree(cls, options):
        '''
        :return: BifurcationTree
        '''
        generationCount = options['Number of generations']
        rootLength = options['Root length']
        rootRadius = 0.5*options['Root diameter']
        forkAngleRadians = radians(options['Fork angle degrees'])
        forkRadiusRatio = options['Fork diameter ratio']
        branchArcRadians = radians(options['Branch arc angle degrees'])
        branchLengthRatio = options['Branch length ratio']
        branchRadiusRatio = options['Branch diameter ratio']
        return BifurcationTree(generationCount, rootLength, rootRadius, forkAngleRadians, forkRadiusRatio, branchArcRadians, branchLengthRatio, branchRadiusRatio)

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base 1D Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        generationCount = options['Number of generations']
        bifurcationTree = cls.generateBifurcationTree(options)
        bifurcationTree.generateZincModel(region)
        return []


class TreeNode:
    '''
    Stores a tree of 1-D bifurcating curves through nested children.
    '''

    def __init__(self, x, d1, r):
        '''
        :param x: coordinates
        :param d1: Parent/primary coordinate derivative
        :param r: Parent/primary radius
        '''
        self._x = x
        self._d1 = [ d1 ]  # list to support child versions
        self._r = [ r ]  # list to support child versions
        self._parent_d1_index = 0
        self._parent_r_index = 0
        self._children = []

    def addChild(self, childTreeNode, d1=None, r=None):
        '''
        :param childTreeNode: Child TreeNode.
        :param d1: Child coordinate derivative, or default None to use parent/primary.
        :param r: Child radius, or default None to use parent/primary.
        '''
        if d1:
            childTreeNode._parent_d1_index = len(self._d1)
            self._d1.append(d1)
        if r:
            childTreeNode._parent_r_index = len(self._r)
            self._r.append(r)
        self._children.append(childTreeNode)

    def getChildCurve(self, childIndex):
        '''
        :param childIndex: Index starting at 0 into child tree node list.
        Return x1, d1, r1, x2, d2, r2
        '''
        assert 0 <= childIndex < len(self._children)
        child = self._children[childIndex]
        return self._x, self._d1[child._parent_d1_index], self._r[child._parent_r_index], child._x, child._d1[0], child._r[0]


class BifurcationTree:
    '''
    Class for generating tree of 1-D bifurcating curves and converting to Zinc model.
    '''

    def __init__(self, generationCount, rootLength, rootRadius, forkAngleRadians, forkRadiusRatio, branchArcRadians, branchLengthRatio, branchRadiusRatio):
        '''
        '''
        self._generationCount = generationCount
        rootDirection = [ 0.0, 0.0, rootLength ]
        self._rootNode = TreeNode([ 0.0, 0.0, 0.0 ], rootDirection, rootRadius)
        self._forkAngleRadians = forkAngleRadians
        self._cosForkAngle = cos(forkAngleRadians)
        self._sinForkAngle = sin(forkAngleRadians)
        self._forkRadiusRatio = forkRadiusRatio
        self._branchArcRadians = branchArcRadians
        self._cosBranchArc = cos(branchArcRadians)
        self._sinBranchArc = sin(branchArcRadians)
        self._branchLengthRatio = branchLengthRatio
        self._branchRadiusRatio = branchRadiusRatio
        self._rootNode.addChild(self._createNodeTree(1, rootDirection, rootDirection, rootRadius*branchRadiusRatio, [ 0.0, 1.0, 0.0 ]))

    def _createNodeTree(self, generation, x1, d1, r, forkNormal):
        '''
        Create node with specified x1, d1, r and recursively add two child nodes until generationCount.
        :param forkNormal: Unit direction normal to d1 and child branches.
        :return: Top node of tree.
        '''
        node = TreeNode(x1, d1, r)
        if generation < self._generationCount:
            branchLength = magnitude(d1)*self._branchLengthRatio
            main = mult(d1, self._cosForkAngle*self._branchLengthRatio)
            side = mult(cross(forkNormal, d1), self._sinForkAngle*self._branchLengthRatio)
            branch1d1 = add(main, side)
            branch2d1 = sub(main, side)
            if self._branchArcRadians > 0.0:
                arcr = branchLength/self._branchArcRadians
                arc2 = mult(branch1d1, arcr/branchLength)
                arc1 = cross(arc2, forkNormal)
                arcc = sub(x1, arc1)
                branch1x2 = add(arcc, add(mult(arc1, self._cosBranchArc), mult(arc2, self._sinBranchArc)))
                branch1d2 = mult(add(mult(arc1, -self._sinBranchArc), mult(arc2, self._cosBranchArc)), branchLength/arcr)
                arc2 = mult(branch2d1, arcr/branchLength)
                arc1 = cross(forkNormal, arc2)
                arcc = sub(x1, arc1)
                branch2x2 = add(arcc, add(mult(arc1, self._cosBranchArc), mult(arc2, self._sinBranchArc)))
                branch2d2 = mult(add(mult(arc1, -self._sinBranchArc), mult(arc2, self._cosBranchArc)), branchLength/arcr)
            else:
                branch1x2 = add(x1, branch1d1)
                branch1d2 = branch1d1
                branch2x2 = add(x1, branch2d1)
                branch2d2 = branch2d1
            branch1Normal = normalize(cross(forkNormal, branch1d2))
            branch2Normal = normalize(cross(forkNormal, branch2d2))
            forkRadius = r*self._forkRadiusRatio
            node.addChild(self._createNodeTree(generation + 1, branch1x2, branch1d2, forkRadius*self._branchRadiusRatio, branch1Normal), branch1d1, forkRadius)
            node.addChild(self._createNodeTree(generation + 1, branch2x2, branch2d2, forkRadius*self._branchRadiusRatio, branch2Normal), branch2d1, forkRadius)
        return node

    def getRootNode(self):
        return self._rootNode

    def generateZincModel(self, region, nextNodeIdentifier=1, nextElementIdentifier=1):
        '''
        Generate Zinc nodes and elements in region to represent tree.
        :return: Final nextNodeIdentifier, nextElementIdentifier.
        '''
        self._fieldmodule = region.getFieldmodule()
        self._nodes = self._fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self._mesh1d = self._fieldmodule.findMeshByDimension(1)
        self._cubicHermiteBasis = self._fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        self._linearBasis = self._fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        self._nodetemplates = {}  # indexed by (d1VersionsCount, rVersionsCount)
        self._elementtemplates = {}  # indexed by start (d1Version, rVersion)
        with ChangeManager(self._fieldmodule):
            self._coordinates = findOrCreateFieldCoordinates(self._fieldmodule)
            self._radius = findOrCreateFieldFiniteElement(self._fieldmodule, "radius", components_count=1, managed=True)
            self._fieldcache = self._fieldmodule.createFieldcache()
            parentNode = None
            self._generateZincModelTree(self._rootNode, parentNode, nextNodeIdentifier, nextElementIdentifier)
        return nextNodeIdentifier, nextElementIdentifier

    def _getZincNodetemplate(self, d1VersionsCount, rVersionsCount):
        '''
        Get node template for specified numbers of d1 and r versions.
        :return: Zinc Nodetemplate
        '''
        templateId = (d1VersionsCount, rVersionsCount)
        nodetemplate = self._nodetemplates.get(templateId)
        if not nodetemplate:
            assert (d1VersionsCount > 0) and (rVersionsCount > 0)
            nodetemplate = self._nodes.createNodetemplate()
            nodetemplate.defineField(self._coordinates)
            nodetemplate.setValueNumberOfVersions(self._coordinates, -1, Node.VALUE_LABEL_D_DS1, d1VersionsCount)
            nodetemplate.defineField(self._radius)
            nodetemplate.setValueNumberOfVersions(self._radius, -1, Node.VALUE_LABEL_VALUE, rVersionsCount)
            self._nodetemplates[templateId] = nodetemplate
        return nodetemplate

    def _getZincElementtemplate(self, d1Version, rVersion):
        '''
        Get node template for specified numbers of d1 and r versions.
        :return: Zinc Elementtemplate
        '''
        templateId = (d1Version, rVersion)
        elementtemplate = self._elementtemplates.get(templateId)
        if not elementtemplate:
            assert (d1Version > 0) and (rVersion > 0)
            elementtemplate = self._mesh1d.createElementtemplate()
            elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
            eftCoordinates = self._mesh1d.createElementfieldtemplate(self._cubicHermiteBasis)
            eftCoordinates.setTermNodeParameter(2, 1, 1, Node.VALUE_LABEL_D_DS1, d1Version)
            elementtemplate.defineField(self._coordinates, -1, eftCoordinates)
            eftRadius = self._mesh1d.createElementfieldtemplate(self._linearBasis)
            eftRadius.setTermNodeParameter(1, 1, 1, Node.VALUE_LABEL_VALUE, rVersion)
            elementtemplate.defineField(self._radius, -1, eftRadius)
            self._elementtemplates[templateId] = elementtemplate
        return elementtemplate

    def _generateZincModelTree(self, treeNode, parentNode, nextNodeIdentifier, nextElementIdentifier):
        '''
        :return: Final nextNodeIdentifier, nextElementIdentifier.
        '''
        d1VersionsCount = len(treeNode._d1)
        rVersionsCount = len(treeNode._r)
        nodetemplate = self._getZincNodetemplate(d1VersionsCount, rVersionsCount)
        node = self._nodes.createNode(nextNodeIdentifier, nodetemplate)
        nextNodeIdentifier += 1
        self._fieldcache.setNode(node)
        self._coordinates.setNodeParameters(self._fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, treeNode._x)
        for i in range(d1VersionsCount):
            self._coordinates.setNodeParameters(self._fieldcache, -1, Node.VALUE_LABEL_D_DS1, i + 1, treeNode._d1[i])
        for i in range(rVersionsCount):
            self._radius.setNodeParameters(self._fieldcache, -1, Node.VALUE_LABEL_VALUE, i + 1, treeNode._r[i])

        if parentNode:
            elementtemplate = self._getZincElementtemplate(treeNode._parent_d1_index + 1, treeNode._parent_r_index + 1)
            element = self._mesh1d.createElement(nextElementIdentifier, elementtemplate)
            nextElementIdentifier += 1
            eftCoordinates = element.getElementfieldtemplate(self._coordinates, -1)
            eftRadius = element.getElementfieldtemplate(self._radius, -1)
            # must set node for both efts
            for eft in ( eftCoordinates, eftRadius ):
                element.setNode(eft, 1, parentNode)
                element.setNode(eft, 2, node)

        for childTreeNode in treeNode._children:
            nextNodeIdentifier, nextElementIdentifier = self._generateZincModelTree(childTreeNode, node, nextNodeIdentifier, nextElementIdentifier)

        return nextNodeIdentifier, nextElementIdentifier
