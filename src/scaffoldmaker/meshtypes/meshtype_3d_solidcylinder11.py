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
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.torso_terms import get_torso_term
from scaffoldmaker.utils import geometry
from scaffoldmaker.utils.shieldmesh import ShieldMesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.interpolation import computeCubicHermiteDerivativeScaling, getCubicHermiteArcLength, interpolateSampleCubicHermite, \
    sampleCubicHermiteCurves, sampleCubicHermiteCurvesSmooth, smoothCubicHermiteDerivativesLine
from opencmiss.utils.zinc.finiteelement import getMaximumElementIdentifier, getMaximumNodeIdentifier
from scaffoldmaker.utils import mirror
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, setEftScaleFactorIds
from scaffoldmaker.utils.cylindermesh import createCylinderMesh3d, CylinderMode

class MeshType_3d_solidcylinder11(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Solid Cylinder 11'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Major radius' : 2.5,
            'Minor radius' : 1.5,
            'Torso height' : 5.0,
            'Number of elements across torso' : 8,
            'Number of elements up torso' : 5,
            'Number of elements along torso' : 5,
            'Number of elements for rim' : 0,
            'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Major radius',
            'Minor radius',
            'Torso height',
            'Number of elements across torso',
            'Number of elements up torso',
            'Number of elements along torso',
            'Number of elements for rim',
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
        majorRadius = options['Major radius']
        minorRadius = options['Minor radius']
        height = options['Torso height']
        elementsCountAcross = options['Number of elements across torso']
        elementsCountUp = options['Number of elements up torso']
        elementsCountRim = options['Number of elements for rim']
        elementsCountAlong = options['Number of elements along torso']
        elementsCountUpRegular = elementsCountUp - 2 - elementsCountRim
        elementsCountAcrossNonRim = elementsCountAcross - 2*elementsCountRim
        elementsCountAround = 2 * elementsCountUpRegular + elementsCountAcrossNonRim
        useCrossDerivatives = options['Use cross derivatives']


        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)


        btGroup = AnnotationGroup(region, get_torso_term("anterior of torso"))
        annotationGroups = [ btGroup ]

        createCylinderMesh3d(fm, coordinates, [0.0, 0.0, 0.0], [0.0, 0.0, height], [majorRadius, 0.0, 0.0], minorRadius,
                             elementsCountAcross, elementsCountUp, elementsCountAlong, 1,
                             1, cylinderMode=CylinderMode.CYLINDER_MODE_FULL, useCrossDerivatives=False)
        createCylinderMesh3d(fm, coordinates, [-3.0, 0.0, 5.0], [-5.0, 0.0, 0.0], [0.0, 0.5, 0.0], 0.5,
                             elementsCountAcross, elementsCountUp, elementsCountAlong, 1,
                             1, cylinderMode=CylinderMode.CYLINDER_MODE_FULL, useCrossDerivatives=False)
        createCylinderMesh3d(fm, coordinates, [3.0, 0.0, 5.0], [5.0, 0.0, 0.0], [0.0, 0.5, 0.0], 0.5,
                             elementsCountAcross, elementsCountUp, elementsCountAlong, 1,
                             1, cylinderMode=CylinderMode.CYLINDER_MODE_FULL, useCrossDerivatives=False)


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
