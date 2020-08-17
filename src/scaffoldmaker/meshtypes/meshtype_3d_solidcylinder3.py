"""
Generates a 3-D cylinder mesh with variable numbers of elements around, along and
across.
"""

from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.cylindermesh2 import CylinderType, CylinderMesh, CylinderMode


class MeshType_3d_solidcylinder3(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Solid Cylinder 3'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Cylinder height' : 5.0,
            'Major radius' : 2.5,
            'Minor radius' : 1.0,
            'Radius reduction rate' : 0.08,
            'Number of elements across' : 8,
            'Number of elements up' : 5,
            'Number of elements along' : 5,
            'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Cylinder height',
            'Major radius',
            'Minor radius',
            'Radius reduction rate',
            'Number of elements across',
            'Number of elements up',
            'Number of elements along',
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
            options['Number of elements through wall'] = 1
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
        height = options['Cylinder height']
        majorRadius = options['Major radius']
        minorRadius = options['Minor radius']
        rate = options['Radius reduction rate']
        elementsCountAcross = options['Number of elements across']
        elementsCountUp = options['Number of elements up']
        elementsCountAlong = options['Number of elements along']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        axis1 = [1.0, 0.0, 0.0]
        axis2 = [0.0, 1.0, 0.0]
        axis3 = [0.0, 0.0, 1.0]

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        cylinder1 = CylinderMesh(fm, coordinates, [0.0, 0.0, 0.0], vector.setMagnitude(axis3, height), vector.setMagnitude(axis1, majorRadius), minorRadius,
                             elementsCountAcross, elementsCountUp, elementsCountAlong,
                             cylinderMode=CylinderMode.CYLINDER_MODE_FULL, cylinderType=CylinderType.CYLIDNER_TRUNCATED_CONE,
                             rate=rate, useCrossDerivatives=False)

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
