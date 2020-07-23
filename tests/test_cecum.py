import unittest
from opencmiss.utils.zinc.finiteelement import evaluateFieldNodesetRange
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.context import Context
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.meshtypes.meshtype_3d_cecum1 import MeshType_3d_cecum1
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace
from testutils import assertAlmostEqualList

class CecumScaffoldTestCase(unittest.TestCase):

    def test_cecum1(self):
        """
        Test creation of cecum scaffold.
        """
        parameterSetNames = MeshType_3d_cecum1.getParameterSetNames()
        self.assertEqual(parameterSetNames, [ "Default", "Pig 1" ])
        options = MeshType_3d_cecum1.getDefaultOptions("Pig 1")
        self.assertEqual(30, len(options))
        self.assertEqual(5, options.get("Number of segments"))
        self.assertEqual(2, options.get("Number of elements around tenia coli"))
        self.assertEqual(8, options.get("Number of elements along segment"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(35.0, options.get("Start inner radius"))
        self.assertEqual(3.0, options.get("Start inner radius derivative"))
        self.assertEqual(38.0, options.get("End inner radius"))
        self.assertEqual(3.0, options.get("End inner radius derivative"))
        self.assertEqual(0.5, options.get("Corner inner radius factor"))
        self.assertEqual(0.25, options.get("Haustrum inner radius factor"))
        self.assertEqual(4.0, options.get("Segment length mid derivative factor"))
        self.assertEqual(3, options.get("Number of tenia coli"))
        self.assertEqual(5.0, options.get("Start tenia coli width"))
        self.assertEqual(0.0, options.get("End tenia coli width derivative"))
        self.assertEqual(2.0, options.get("Wall thickness"))
        ostiumOptions = options['Ileocecal junction']
        ostiumSettings = ostiumOptions.getScaffoldSettings()
        self.assertEqual(1, ostiumSettings.get("Number of vessels"))
        self.assertEqual(8, ostiumSettings.get("Number of elements around ostium"))
        self.assertEqual(1, ostiumSettings.get("Number of elements through wall"))
        self.assertEqual(20.0, ostiumSettings.get("Ostium diameter"))
        self.assertEqual(10.0, ostiumSettings.get("Vessel inner diameter"))
        self.assertEqual(60, options.get("Ileocecal junction angular position degrees"))
        self.assertEqual(0.5, options.get("Ileocecal junction position along factor"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = MeshType_3d_cecum1.generateBaseMesh(region, options)
        self.assertEqual(2, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(1492, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(5617, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(6767, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(2642, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-49.01658984455258, -46.89686037622053, -2.343256155753525], 1.0E-6)
        assertAlmostEqualList(self, maximums, [42.18085849205387, 54.89264119402881, 180.0], 1.0E-6)

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
        self.assertAlmostEqual(surfaceArea, 65960.02821062482, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 127893.74708048582, delta=1.0E-6)

if __name__ == "__main__":
    unittest.main()
