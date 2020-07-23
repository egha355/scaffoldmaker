import copy
import unittest
from opencmiss.utils.zinc.finiteelement import evaluateFieldNodesetRange
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.context import Context
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_smallintestine1 import MeshType_3d_smallintestine1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace, exnodeStringFromNodeValues
from testutils import assertAlmostEqualList

class SmallIntestineScaffoldTestCase(unittest.TestCase):

    def test_smallintestine1(self):
        """
        Test creation of small intestine scaffold.
        """
        parameterSetNames = MeshType_3d_smallintestine1.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Mouse 1"])
        centralPathDefaultScaffoldPackages = {
            'Test line': ScaffoldPackage(MeshType_1d_path1, {
                'scaffoldSettings': {
                    'Coordinate dimensions': 3,
                    'Length': 1.0,
                    'Number of elements': 3
                },
                'meshEdits': exnodeStringFromNodeValues(
                    [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                     Node.VALUE_LABEL_D2_DS1DS2], [
                        [[-2.3, 18.5, -4.4], [-4.2, -0.8, 3.7], [0.0, 5.0, 0.0], [0.0, 0.0, 0.5]],
                        [[-8.6, 16.3, -0.4], [-7.1, -2.7, 1.6], [0.0, 5.0, 0.0], [0.0, 0.0, 0.5]],
                        [[-18.3, 12.6, -1.5], [-6.4, -1.7, -3.8], [0.0, 5.0, 0.0], [0.0, 0.0, 0.5]],
                        [[-15.6, 13.7, -6.1], [7.0, 2.1, -1.8], [0.0, 5.0, 0.0], [0.0, 0.0, 0.5]]])
            })
        }
        centralPathOption = centralPathDefaultScaffoldPackages['Test line']
        options = MeshType_3d_smallintestine1.getDefaultOptions("Mouse 1")
        options['Central path'] = copy.deepcopy(centralPathOption)
        options['Number of segments'] = 4
        options['Duodenum length'] = 5.0
        options['Jejunum length'] = 15.0
        options['Ileum length'] = 5.0
        self.assertEqual(19, len(options))
        centralPath = options['Central path']
        self.assertEqual(4, options.get("Number of segments"))
        self.assertEqual(8, options.get("Number of elements around"))
        self.assertEqual(4, options.get("Number of elements along segment"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(5.0, options.get("Duodenum length"))
        self.assertEqual(5.0, options.get("Ileum length"))
        self.assertEqual(0.6, options.get("Duodenum inner radius"))
        self.assertEqual(1.0, options.get("Jejunum-ileum inner radius"))
        self.assertEqual(0.1, options.get("Wall thickness"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())

        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx = extractPathParametersFromRegion(tmpRegion)[0]
        self.assertEqual(4, len(cx))
        assertAlmostEqualList(self, cx[0], [  -2.3, 18.5, -4.4 ], 1.0E-6)
        assertAlmostEqualList(self, cx[2], [ -18.3, 12.6, -1.5 ], 1.0E-6)
        del tmpRegion

        annotationGroups = MeshType_3d_smallintestine1.generateBaseMesh(region, options)
        self.assertEqual(4, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(128, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(520, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(664, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(272, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [ -20.06978981419564, 11.406595205949705, -7.1653294859433965 ], 1.0E-6)
        assertAlmostEqualList(self, maximums, [ -1.8300388314851923, 19.193885338090105, 0.9772071374844936 ], 1.0E-6)

        flatCoordinates = fieldmodule.findFieldByName("flat coordinates").castFiniteElement()
        self.assertTrue(flatCoordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(flatCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [ -1.39038154442654, 0.0, 0.0 ], 1.0E-6)
        assertAlmostEqualList(self, maximums, [ 4.891237158967401, 25.293706698841913, 0.1 ], 1.0E-6)

        textureCoordinates = fieldmodule.findFieldByName("texture coordinates").castFiniteElement()
        minimums, maximums = evaluateFieldNodesetRange(textureCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [ 0.0, 0.0, 0.0 ], 1.0E-6)
        assertAlmostEqualList(self, maximums, [ 0.875, 1.0, 1.0 ], 1.0E-6)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            faceMeshGroup = createFaceMeshGroupExteriorOnFace(fieldmodule, Element.FACE_TYPE_XI3_1)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, faceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
            flatSurfaceAreaField = fieldmodule.createFieldMeshIntegral(one, flatCoordinates, faceMeshGroup)
            flatSurfaceAreaField.setNumbersOfPoints(4)
            textureSurfaceAreaField = fieldmodule.createFieldMeshIntegral(one, textureCoordinates, faceMeshGroup)
            textureSurfaceAreaField.setNumbersOfPoints(4)
            textureVolumeField = fieldmodule.createFieldMeshIntegral(one, textureCoordinates, mesh3d)
            textureVolumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 171.27464080337143, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 16.35219225882822, delta=1.0E-6)
        result, flatSurfaceArea = flatSurfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(flatSurfaceArea, 171.37026123844635, delta=1.0E-3)
        result, textureSurfaceArea = textureSurfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(textureSurfaceArea, 1.0, delta=1.0E-6)
        result, textureVolume = textureVolumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(textureVolume, 1.0, delta=1.0E-6)

if __name__ == "__main__":
    unittest.main()
