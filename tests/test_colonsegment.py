import unittest
from opencmiss.utils.zinc.finiteelement import evaluateFieldNodesetRange
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.context import Context
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import MeshType_3d_colonsegment1
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace
from testutils import assertAlmostEqualList

class ColonSegmentScaffoldTestCase(unittest.TestCase):

    def test_humancolonsegment1(self):
        """
        Test creation of human colon segment scaffold.
        """
        parameterSetNames = MeshType_3d_colonsegment1.getParameterSetNames()
        self.assertEqual(parameterSetNames, [ "Default", "Human 1", "Mouse 1", "Pig 1" ])
        options = MeshType_3d_colonsegment1.getDefaultOptions("Human 1")
        self.assertEqual(27, len(options))
        self.assertEqual(0.0, options.get("Start phase"))
        self.assertEqual(2, options.get("Number of elements around tenia coli"))
        self.assertEqual(4, options.get("Number of elements along segment"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(43.5, options.get("Start inner radius"))
        self.assertEqual(0.0, options.get("Start inner radius derivative"))
        self.assertEqual(33.0, options.get("End inner radius"))
        self.assertEqual(0.0, options.get("End inner radius derivative"))
        self.assertEqual(3.0, options.get("Segment length mid derivative factor"))
        self.assertEqual(50.0, options.get("Segment length"))
        self.assertEqual(3, options.get("Number of tenia coli"))
        self.assertEqual(10.0, options.get("Start tenia coli width"))
        self.assertEqual(0.0, options.get("End tenia coli width derivative"))
        self.assertEqual(1.6, options.get("Wall thickness"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = MeshType_3d_colonsegment1.generateBaseMesh(region, options)
        self.assertEqual(3, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(144, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(576, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(747, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(315, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [ -2.172286248499807e-15, -58.95670186936737, -55.54662267827035 ], 1.0E-6)
        assertAlmostEqualList(self, maximums, [ 50.0, 50.52621132610023, 55.54662267827035 ], 1.0E-6)

        flatCoordinates = fieldmodule.findFieldByName("flat coordinates").castFiniteElement()
        self.assertTrue(flatCoordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(flatCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [ 0.0, 0.0, 0.0 ], 1.0E-6)
        assertAlmostEqualList(self, maximums, [397.2736607240895, 50.0, 3.2000000000000006], 1.0E-6)

        textureCoordinates = fieldmodule.findFieldByName("texture coordinates").castFiniteElement()
        minimums, maximums = evaluateFieldNodesetRange(textureCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [ 0.0, 0.0, 0.0 ], 1.0E-6)
        assertAlmostEqualList(self, maximums, [ 0.9887737064410717, 1.0, 2.0 ], 1.0E-6)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            faceMeshGroup = createFaceMeshGroupExteriorOnFace(fieldmodule, Element.FACE_TYPE_XI3_1)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, faceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 21129.564192298032, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 40783.41115796618, delta=1.0E-6)

    def test_mousecolonsegment1(self):
        """
        Test creation of mouse colon segment scaffold.
        """
        options = MeshType_3d_colonsegment1.getDefaultOptions("Mouse 1")
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = MeshType_3d_colonsegment1.generateBaseMesh(region, options)
        self.assertEqual(2, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(40, mesh3d.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        flatCoordinates = fieldmodule.findFieldByName("flat coordinates").castFiniteElement()
        self.assertTrue(flatCoordinates.isValid())
        textureCoordinates = fieldmodule.findFieldByName("texture coordinates").castFiniteElement()
        self.assertTrue(textureCoordinates.isValid())

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            faceMeshGroup = createFaceMeshGroupExteriorOnFace(fieldmodule, Element.FACE_TYPE_XI3_1)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, faceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            flatSurfaceAreaField = fieldmodule.createFieldMeshIntegral(one, flatCoordinates, faceMeshGroup)
            flatSurfaceAreaField.setNumbersOfPoints(4)
            textureSurfaceAreaField = fieldmodule.createFieldMeshIntegral(one, textureCoordinates, faceMeshGroup)
            textureSurfaceAreaField.setNumbersOfPoints(4)
            textureVolumeField = fieldmodule.createFieldMeshIntegral(one, textureCoordinates, mesh3d)
            textureVolumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 468.17062489996886, delta=1.0E-6)
        result, flatSurfaceArea = flatSurfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(flatSurfaceArea, surfaceArea, delta=1.0E-3)
        result, textureSurfaceArea = textureSurfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(textureSurfaceArea, 1.0, delta=1.0E-6)
        result, textureVolume = textureVolumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(textureVolume, 1.0, delta=1.0E-6)

if __name__ == "__main__":
    unittest.main()
