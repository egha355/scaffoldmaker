"""
Generates a 3-D unit tube mesh with variable numbers of elements around, along and
through wall, plus variable wall thickness for unit diameter.
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


class MeshType_3d_wholebody1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Whole Body 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements around' : 4,
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
        elementsCountAround = options['Number of elements around']
        elementsCountAlong = options['Number of elements along']
        elementsCountThroughWall = options['Number of elements through wall']
        useCrossDerivatives = options['Use cross derivatives']

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
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        mesh = fm.findMeshByDimension(3)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft2 = tricubichermite.createEftBasic()
        eft = tricubichermite.createEftWedgeRadial(0 * 100, 1* 100)
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        elementtemplate2 = mesh.createElementtemplate()
        elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate2.defineField(coordinates, -1, eft2)

        cache = fm.createFieldcache()

        # create nodes
        nodeIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.5/elementsCountThroughWall , 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 1.0 / elementsCountAlong ]
        dx_ds3 = [ 0.0, 0.5/elementsCountThroughWall, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        wallThicknessPerElement = 0.5/elementsCountThroughWall

        for n2 in range(elementsCountAlong+1):
            x = [0.0, 0.0, n2 * 1.0 / elementsCountAlong]
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
            nodeIdentifier += 1

        for n3 in range(elementsCountThroughWall):
            radius = 0.5*((n3+1)/elementsCountThroughWall)
            for n2 in range(elementsCountAlong + 1):
                x[2] = n2 / elementsCountAlong
                for n1 in range(elementsCountAround):
                    radiansAround = n1*radiansPerElementAround
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    x[0] = radius*cosRadiansAround
                    x[1] = radius*sinRadiansAround
                    dx_ds1[0] = radiansPerElementAround*radius*-sinRadiansAround
                    dx_ds1[1] = radiansPerElementAround*radius*cosRadiansAround
                    dx_ds3[0] = wallThicknessPerElement*cosRadiansAround
                    dx_ds3[1] = wallThicknessPerElement*sinRadiansAround
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    if useCrossDerivatives:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                    nodeIdentifier = nodeIdentifier + 1

        # create elements
        elementIdentifier = 1
        now = (elementsCountAlong + 1)*elementsCountAround
        for e3 in range(elementsCountThroughWall):
            for e2 in range(elementsCountAlong):
                for e1 in range(elementsCountAround):
                    if e3 == 0:
                        element = mesh.createElement(elementIdentifier, elementtemplate)
                        bni1 = e2+1
                        bni2 = elementsCountAlong+1+e1+1+e2*elementsCountAround
                        bni3 = elementsCountAlong+1+(e1+1) % elementsCountAround+1+e2*elementsCountAround
                        nodeIdentifiers = [ bni1, bni1+1, bni2, bni3, bni2+elementsCountAround, bni3+elementsCountAround ]
                        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                        elementIdentifier = elementIdentifier + 1
                    else:
                        element = mesh.createElement(elementIdentifier, elementtemplate2)
                        bni11 = e3*now + e2*elementsCountAround + e1 + 1
                        bni12 = e3*now + e2*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                        bni21 = e3*now + (e2 + 1)*elementsCountAround + e1 + 1
                        bni22 = e3*now + (e2 + 1)*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                        # nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + now, bni12 + now, bni21 + now, bni22 + now ]
                        bni1 = elementsCountAlong+1+1+e1+e2*elementsCountAround+(e3-1)*now
                        bni2 = elementsCountAlong+1+1+(e1+1)%elementsCountAround+e2*elementsCountAround+(e3-1)*now
                        nodeIdentifiers = [ bni1, bni2, bni1+elementsCountAround, bni2+elementsCountAround, bni1+now, bni2+now,bni1+elementsCountAround+now, bni2+elementsCountAround+now]
                        result = element.setNodesByIdentifier(eft2, nodeIdentifiers)
                        elementIdentifier = elementIdentifier + 1
        fm.endChange()

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
