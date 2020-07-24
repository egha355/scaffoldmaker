import unittest
import copy
from opencmiss.utils.zinc.finiteelement import evaluateFieldNodesetRange, findNodeWithName
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.context import Context
from opencmiss.zinc.field import Field
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, getAnnotationGroupForTerm
from scaffoldmaker.annotation.heart_terms import get_heart_term
from scaffoldmaker.meshtypes.meshtype_3d_heart1 import MeshType_3d_heart1
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from testutils import assertAlmostEqualList

class HeartScaffoldTestCase(unittest.TestCase):

    def test_heart1(self):
        """
        Test creation of heart scaffold.
        """
        scaffold = MeshType_3d_heart1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, [ "Default", "Human 1", "Mouse 1", "Pig 1", "Rat 1",
            "Unit Human 1", "Unit Mouse 1", "Unit Pig 1", "Unit Rat 1" ]);
        options = scaffold.getDefaultOptions("Human 1")
        self.assertEqual(117, len(options))
        self.assertEqual(0.9, options.get("LV outer height"))
        self.assertEqual(80.0, options.get("Unit scale"))
        self.assertEqual(7, options.get("Number of elements around LV free wall"))
        self.assertEqual(7, options.get("Number of elements around RV free wall"))
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateMesh(region, options)
        self.assertEqual(27, len(annotationGroups))
        fieldmodule = region.getFieldmodule()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(289, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(1117, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(1356, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(528, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check coordinates range, epicardium surface area and volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [ -50.7876375290527, -57.76496144495844, -91.6 ], 1.0E-6)
        assertAlmostEqualList(self, maximums, [ 43.810947610743156, 39.03925080604259, 40.717553621061704 ], 1.0E-6)
        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            epicardiumGroup = fieldmodule.findFieldByName('epicardium').castGroup()
            self.assertTrue(epicardiumGroup.isValid())
            epicardiumMeshGroup = epicardiumGroup.getFieldElementGroup(mesh2d).getMeshGroup()
            self.assertTrue(epicardiumMeshGroup.isValid())
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, epicardiumMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 36545.32654122452, delta=1.0E-3)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 222019.91858237237, delta=1.0E-3)

        # check some annotationGroups:
        expectedSizes3d = {
            "left ventricle myocardium" : 94,
            "right ventricle myocardium" : 79,
            "interventricular septum" : 30,
            "left atrium myocardium" : 88,
            "right atrium myocardium" : 62,
            "interatrial septum" : 17
            }
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_heart_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)
        expectedSizes2d = {
            "endocardium of left ventricle" : 74,
            "endocardium of right ventricle" : 59,
            "endocardium of left atrium" : 82,
            "endocardium of right atrium" : 56,
            "epicardium" : 197
            }
        for name in expectedSizes2d:
            group = getAnnotationGroupForTerm(annotationGroups, get_heart_term(name))
            size = group.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name], size, name)

        # refine 2x2x2 and check result
        # first remove any surface annotation groups as they are re-added by defineFaceAnnotations
        removeAnnotationGroups = []
        for annotationGroup in annotationGroups:
            if annotationGroup.getMeshGroup(mesh3d).getSize() == 0:
                removeAnnotationGroups.append(annotationGroup)
        for annotationGroup in removeAnnotationGroups:
            annotationGroups.remove(annotationGroup)
        self.assertEqual(18, len(annotationGroups))

        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine number of elements surface'] = 2
        options['Refine number of elements through LV wall'] = 2
        options['Refine number of elements through wall'] = 2
        meshrefinement = MeshRefinement(region, refineRegion, annotationGroups)
        scaffold.refineMesh(meshrefinement, options)
        annotationGroups = meshrefinement.getAnnotationGroups()

        refineFieldmodule.defineAllFaces()
        oldAnnotationGroups = copy.copy(annotationGroups)
        for annotationGroup in annotationGroups:
            annotationGroup.addSubelements()
        scaffold.defineFaceAnnotations(refineRegion, options, annotationGroups)
        for annotation in annotationGroups:
            if annotation not in oldAnnotationGroups:
                annotationGroup.addSubelements()
        self.assertEqual(27, len(annotationGroups))

        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(2236, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(7676, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(8623, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(3183, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_heart_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name]*8, size, name)
        for name in expectedSizes2d:
            group = getAnnotationGroupForTerm(annotationGroups, get_heart_term(name))
            size = group.getMeshGroup(mesh2d).getSize()
            if name == "epicardium":
                # fibrous ring is only refined by 1 through its thickness
                fibrousRingReduction = (options['Refine number of elements surface'] - 1)*options['Refine number of elements surface']*(
                    options['Number of elements around left atrium free wall'] + options['Number of elements around right atrium free wall'])
                self.assertEqual(expectedSizes2d[name]*4 - fibrousRingReduction, size, name)
            else:
                self.assertEqual(expectedSizes2d[name]*4, size, name)

        # test finding a marker in refined scaffold
        markerGroup = refineFieldmodule.findFieldByName("marker").castGroup()
        refinedNodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNodes = markerGroup.getFieldNodeGroup(refinedNodes).getNodesetGroup()
        self.assertEqual(5, markerNodes.getSize())
        markerName = refineFieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = refineFieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        cache = refineFieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "APEX")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertEqual(5, element.getIdentifier())
        assertAlmostEqualList(self, xi, [ 0.0, 0.0, 1.0 ], 1.0E-10)

if __name__ == "__main__":
    unittest.main()
