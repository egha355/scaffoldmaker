"""
Generates 3-D mesh of whole body.
"""

from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.torso_terms import get_torso_term
from scaffoldmaker.utils.cylindermesh2 import CylinderMesh, CylinderMode, CylinderType
from scaffoldmaker.utils import vector

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
            'Torso major radius' : 2.5,
            'Torso minor radius' : 1.5,
            'Torso height' : 5.0,
            'Torso radius reduction rate' : 0.08,
            'Number of elements across torso' : 8,
            'Number of elements up torso' : 5,
            'Number of elements along torso' : 5,
            'Number of elements for rim' : 0,
            'Arm radius' : 0.5,
            'Arm length' : 5.0,
            'Arm centre height on axis3' : 5.0,
            'Arm base distance from centre on axis3' : 3.0,
            'Arm angle degrees': 90,
            'Leg radius': 0.75,
            'Leg length': 7.5,
            'Leg centre height on axis3': 0.5,
            'Leg base distance from centre on axis3': 2.0,
            'Leg angle degrees': 30,
            'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Torso major radius',
            'Torso minor radius',
            'Torso height',
            'Torso radius reduction rate',
            'Number of elements across torso',
            'Number of elements up torso',
            'Number of elements along torso',
            'Number of elements for rim',
            'Arm radius',
            'Arm length',
            'Arm centre height on axis3',
            'Arm base distance from centre on axis3',
            'Arm angle degrees',
            'Leg radius',
            'Leg length',
            'Leg centre height on axis3',
            'Leg base distance from centre on axis3',
            'Leg angle degrees',
            'Use cross derivatives',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements across torso',
            'Number of elements along torso']:
            if options[key] < 1:
                options[key] = 1
        if (options['Number of elements across torso'] < 2) :
            options['Number of elements across torso'] = 8


    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        majorRadius = options['Torso major radius']
        minorRadius = options['Torso minor radius']
        height = options['Torso height']
        rate = options['Torso radius reduction rate']
        elementsCountAcross = options['Number of elements across torso']
        elementsCountUp = options['Number of elements up torso']
        elementsCountRim = options['Number of elements for rim']
        elementsCountAlong = options['Number of elements along torso']
        elementsCountUpRegular = elementsCountUp - 2 - elementsCountRim
        elementsCountAcrossNonRim = elementsCountAcross - 2*elementsCountRim
        elementsCountAround = 2 * elementsCountUpRegular + elementsCountAcrossNonRim
        useCrossDerivatives = options['Use cross derivatives']
        armRadius = options['Arm radius']
        armLength = options['Arm length']
        armh = options['Arm centre height on axis3']
        armd = options['Arm base distance from centre on axis3']
        armt = math.radians(options['Arm angle degrees'])
        legRadius = options['Leg radius']
        legLength = options['Leg length']
        legh = options['Leg centre height on axis3']
        legd = options['Leg base distance from centre on axis3']
        legt = math.radians(options['Leg angle degrees'])

        axis1 = [1.0, 0.0, 0.0]
        axis2 = [0.0, 1.0, 0.0]
        axis3 = [0.0, 0.0, 1.0]

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)


        btGroup = AnnotationGroup(region, get_torso_term("anterior of torso"))
        annotationGroups = [ btGroup ]

        # set arms cylinder properties
        larmBaseCentre = vector.addVectors(axis1, axis3, armd*math.sin(armt), armh-armd*math.cos(armt))
        rarmBaseCentre = vector.addVectors(axis1, axis3, armd*math.sin(-armt), armh-armd*math.cos(-armt))
        larmAlongAxis  = vector.addVectors(axis1, axis3, armLength*math.sin(armt), -armLength*math.cos(armt))
        rarmAlongAxis  = vector.addVectors(axis1, axis3, armLength*math.sin(-armt), -armLength*math.cos(-armt))
        larmMajorAxis   = vector.addVectors(axis1, axis3, armRadius*math.cos(armt), armRadius*math.sin(armt))
        rarmMajorAxis   = vector.addVectors(axis1, axis3, armRadius*math.cos(-armt), armRadius*math.sin(-armt))
        # set legs cylinders properties
        llegBaseCentre = vector.addVectors(axis1, axis3, legd*math.sin(legt), legh-legd*math.cos(legt))
        rlegBaseCentre = vector.addVectors(axis1, axis3, legd*math.sin(-legt), legh-legd*math.cos(-legt))
        llegAlongAxis  = vector.addVectors(axis1, axis3, legLength*math.sin(legt), -legLength*math.cos(legt))
        rlegAlongAxis  = vector.addVectors(axis1, axis3, legLength*math.sin(-legt), -legLength*math.cos(-legt))
        llegMajorAxis   = vector.addVectors(axis1, axis3, legRadius*math.cos(legt), legRadius*math.sin(legt))
        rlegMajorAxis   = vector.addVectors(axis1, axis3, legRadius*math.cos(-legt), legRadius*math.sin(-legt))

        torso = CylinderMesh(fm, coordinates, [0.0, 0.0, 0.0], vector.setMagnitude(axis3, height), vector.setMagnitude(axis1, majorRadius), minorRadius,
                             elementsCountAcross, elementsCountUp, elementsCountAlong,
                             cylinderMode=CylinderMode.CYLINDER_MODE_FULL, cylinderType=CylinderType.CYLIDNER_TRUNCATED_CONE,
                             rate=rate, useCrossDerivatives=False)
        larm = CylinderMesh(fm, coordinates, larmBaseCentre, larmAlongAxis, larmMajorAxis, armRadius,
                             elementsCountAcross, elementsCountUp, elementsCountAlong,
                             cylinderMode=CylinderMode.CYLINDER_MODE_FULL, useCrossDerivatives=False)
        rarm = CylinderMesh(fm, coordinates, rarmBaseCentre, rarmAlongAxis, rarmMajorAxis, armRadius,
                             elementsCountAcross, elementsCountUp, elementsCountAlong,
                             cylinderMode=CylinderMode.CYLINDER_MODE_FULL, useCrossDerivatives=False)
        lleg = CylinderMesh(fm, coordinates, llegBaseCentre, llegAlongAxis, llegMajorAxis, legRadius,
                             elementsCountAcross, elementsCountUp, elementsCountAlong,
                             cylinderMode=CylinderMode.CYLINDER_MODE_FULL, useCrossDerivatives=False)
        rleg = CylinderMesh(fm, coordinates, rlegBaseCentre, rlegAlongAxis, rlegMajorAxis, legRadius,
                             elementsCountAcross, elementsCountUp, elementsCountAlong,
                             cylinderMode=CylinderMode.CYLINDER_MODE_FULL, useCrossDerivatives=False)



        fm.endChange()
        return annotationGroups

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
        anteriorGroup = getAnnotationGroupForTerm(annotationGroups, get_torso_term("anterior of torso"))
        posteriorGroup = getAnnotationGroupForTerm(annotationGroups, get_torso_term("posterior of torso"))

        mesh2d = fm.findMeshByDimension(2)
        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))

        is_anterior = anteriorGroup.getFieldElementGroup(mesh2d)
        is_anterior_faces = fm.createFieldAnd(is_anterior, is_exterior_face_xi3_1)


        is_posterior = posteriorGroup.getFieldElementGroup(mesh2d)
        is_posterior_faces = fm.createFieldAnd(is_posterior, is_exterior_face_xi3_1)

        anteriorFaces = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_torso_term("anterior of torso"))
        anteriorFaces.getMeshGroup(mesh2d).addElementsConditional(is_anterior_faces)
        posteriorFaces = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_torso_term("posterior of torso"))
        posteriorFaces.getMeshGroup(mesh2d).addElementsConditional(is_posterior_faces)
