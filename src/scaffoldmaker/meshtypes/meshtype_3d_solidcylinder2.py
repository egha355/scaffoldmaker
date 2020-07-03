"""
Generates a 3-D unit cylinder mesh with variable numbers of elements around, along and
through wall.
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


class MeshType_3d_solidcylinder2(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Solid Cylinder 2'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements around' : 8,
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
        if (options['Number of elements through wall'] < 2) :
            options['Number of elements around'] = 2
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
        eft = tricubichermite.createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        cache = fm.createFieldcache()

        # create nodes
        nodeIdentifier = 1
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0 , 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 1.0 / elementsCountAlong ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]

        dist = 0.2
        radius = 0.5
        for nz in range(elementsCountAlong+1):
            x = [0.0, 0.0, nz * elementsCountAlong]
            dx_ds1 = [dist, 0.0, 0.0]
            dx_ds3 = [0.0, dist, 0.0]
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
            nodeIdentifier += 1

        # inner square
        cx = [1, 1, 0, -1, -1, -1, 0, 1]
        cy = [0, 1, 1, 1, 0, -1, -1, -1]
        for nz in range(elementsCountAlong+1):
            for nr in range(8):
                x = [dist*cx[nr], dist*cy[nr], nz * elementsCountAlong]
                dx_ds1 = [dist, 0.0, 0.0]
                dx_ds3 = [0.0, dist, 0.0]
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1

        # outer circle
        radiansPerElementAround = 2.0 * math.pi / elementsCountAround
        for nz in range(elementsCountAlong+1):
            x[2] = nz * 1.0/elementsCountAlong
            for nr in range(elementsCountThroughWall):
                for nt in range(elementsCountAround):
                    radiansAround = nt * radiansPerElementAround
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    x[0] = radius*cosRadiansAround
                    x[1] = radius*sinRadiansAround
                    dx_ds3 = [-radius*radiansPerElementAround*sinRadiansAround, radius*radiansPerElementAround*cosRadiansAround, 0.0]
                    dx_ds1 = [(radius-dist)*cosRadiansAround, (radius-dist)*sinRadiansAround, 0.0]
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    nodeIdentifier += 1


        # create elements
        # inner cubic elements
        elementIdentifier = 1
        nodeIdentifiers = [1, 3, 2, 11, 5, 4, 13, 12]
        element = mesh.createElement(elementIdentifier, elementtemplate)
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1
        nodeIdentifiers = [7, 1, 15, 2, 6, 5, 14, 13]
        element = mesh.createElement(elementIdentifier, elementtemplate)
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1
        nodeIdentifiers = [8, 9, 16, 17, 7, 1, 15, 2]
        element = mesh.createElement(elementIdentifier, elementtemplate)
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1
        nodeIdentifiers = [9, 10, 17, 18, 1, 3, 2, 11]
        element = mesh.createElement(elementIdentifier, elementtemplate)
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        # outer elements
        # now = (elementsCountAlong + 1)*elementsCountAround
        # for e3 in range(elementsCountThroughWall+1):
        for e2 in range(elementsCountAlong):
            for e1 in range(elementsCountAround):
                bni1 = elementsCountAlong+1+e1+1
                bni2 = elementsCountAlong+1+(e1+1)%elementsCountAround+1
                # nodeIdentifiers = [bni1, bni2, bni1+elementsCountAround, bni2+elementsCountAround, bni1+2*elementsCountAround, bni2+2*elementsCountAround, bni1+3*elementsCountAround, bni2+3*elementsCountAround]
                nodeIdentifiers = [bni1, bni1+2*elementsCountAround, bni1+elementsCountAround, bni1+3*elementsCountAround,bni2 , bni2+2*elementsCountAround,  bni2+elementsCountAround, bni2+3*elementsCountAround]
                element = mesh.createElement(elementIdentifier, elementtemplate)
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

        #             if e3 == 0:
        #                 element = mesh.createElement(elementIdentifier, elementtemplate)
        #                 bni1 = e2+1
        #                 bni2 = elementsCountAlong+1+e1+1+e2*elementsCountAround
        #                 bni3 = elementsCountAlong+1+(e1+1) % elementsCountAround+1+e2*elementsCountAround
        #                 nodeIdentifiers = [ bni1, bni1+1, bni2, bni3, bni2+elementsCountAround, bni3+elementsCountAround ]
        #                 result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        #                 elementIdentifier = elementIdentifier + 1
        #             else:
        #                 element = mesh.createElement(elementIdentifier, elementtemplate2)
        #                 bni11 = e3*now + e2*elementsCountAround + e1 + 1
        #                 bni12 = e3*now + e2*elementsCountAround + (e1 + 1)%elementsCountAround + 1
        #                 bni21 = e3*now + (e2 + 1)*elementsCountAround + e1 + 1
        #                 bni22 = e3*now + (e2 + 1)*elementsCountAround + (e1 + 1)%elementsCountAround + 1
        #                 # nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + now, bni12 + now, bni21 + now, bni22 + now ]
        #                 bni1 = elementsCountAlong+1+1+e1+e2*elementsCountAround+(e3-1)*now
        #                 bni2 = elementsCountAlong+1+1+(e1+1)%elementsCountAround+e2*elementsCountAround+(e3-1)*now
        #                 nodeIdentifiers = [ bni1, bni2, bni1+elementsCountAround, bni2+elementsCountAround, bni1+now, bni2+now,bni1+elementsCountAround+now, bni2+elementsCountAround+now]
        #                 result = element.setNodesByIdentifier(eft2, nodeIdentifiers)
        #                 elementIdentifier = elementIdentifier + 1
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
