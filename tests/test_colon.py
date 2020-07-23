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
from scaffoldmaker.meshtypes.meshtype_3d_colon1 import MeshType_3d_colon1
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import MeshType_3d_colonsegment1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace, exnodeStringFromNodeValues
from testutils import assertAlmostEqualList

class ColonScaffoldTestCase(unittest.TestCase):

    def test_colon1(self):
        """
        Test creation of colon scaffold.
        """
        parameterSetNames = MeshType_3d_colon1.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Human 1", "Human 2", "Mouse 1", "Mouse 2", "Pig 1", "Pig 2"])
        centralPathDefaultScaffoldPackages = {
            'Test line': ScaffoldPackage(MeshType_1d_path1, {
                'scaffoldSettings': {
                    'Coordinate dimensions': 3,
                    'Length': 1.0,
                    'Number of elements': 1
                },
                'meshEdits': exnodeStringFromNodeValues(
                    [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                     Node.VALUE_LABEL_D2_DS1DS2], [
                        [[163.7, -25.2, 12.2], [-21.7, 50.1, -18.1], [0.0, 0.0, 5.0], [0.0, 0.0, 0.5]],
                        [[117.2, 32.8, -2.6], [-64.3, 34.4, -3.9], [0.0, 0.0, 5.0], [0.0, 0.0, 0.5]]])
            })
        }
        centralPathOption = centralPathDefaultScaffoldPackages['Test line']
        segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName='Human 1')
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Segment profile': segmentProfileOption,
            'Number of segments': 3,
            'Start phase': 0.0,
            'Proximal length': 25.0,
            'Transverse length': 25.0,
            'Distal length': 25.0,
            'Proximal inner radius': 20.0,
            'Proximal tenia coli width': 8.0,
            'Proximal-transverse inner radius': 18.0,
            'Proximal-transverse tenia coli width': 6.0,
            'Transverse-distal inner radius': 16.0,
            'Transverse-distal tenia coli width': 5.0,
            'Distal inner radius': 15.0,
            'Distal tenia coli width': 5.0,
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements along': 1,
            'Refine number of elements through wall': 1
        }
        self.assertEqual(19, len(options))
        centralPath = options['Central path']
        segmentProfile = options.get("Segment profile")
        segmentSettings = segmentProfile.getScaffoldSettings()
        self.assertEqual(8, segmentSettings.get("Number of elements around haustrum"))
        self.assertEqual(0.5, segmentSettings.get("Corner inner radius factor"))
        self.assertEqual(0.5, segmentSettings.get("Haustrum inner radius factor"))
        self.assertEqual(0.5, segmentSettings.get("Segment length end derivative factor"))
        self.assertEqual(3, segmentSettings.get("Number of tenia coli"))
        self.assertEqual(1.6, segmentSettings.get("Tenia coli thickness"))
        self.assertEqual(3, options.get("Number of segments"))
        self.assertEqual(0.0, options.get("Start phase"))
        self.assertEqual(25.0, options.get("Transverse length"))
        self.assertEqual(20.0, options.get("Proximal inner radius"))
        self.assertEqual(6.0, options.get("Proximal-transverse tenia coli width"))
        self.assertEqual(16.0, options.get("Transverse-distal inner radius"))
        self.assertEqual(5.0, options.get("Distal tenia coli width"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())

        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx = extractPathParametersFromRegion(tmpRegion)[0]
        self.assertEqual(2, len(cx))
        assertAlmostEqualList(self, cx[0], [ 163.7, -25.2, 12.2 ], 1.0E-6)
        assertAlmostEqualList(self, cx[1], [ 117.2, 32.8, -2.6 ], 1.0E-6)
        del tmpRegion

        annotationGroups = MeshType_3d_colon1.generateBaseMesh(region, options)
        self.assertEqual(7, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        if annotationGroups is not None:
            for annotationGroup in annotationGroups:
                annotationGroup.addSubelements()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(432, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(1656, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(2043, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(819, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [ 108.0250647983898, -36.876103983560014,  -25.89741158325083 ], 1.0E-6)
        assertAlmostEqualList(self, maximums, [ 185.46457506220003, 48.101157490744, 34.995316052158934 ], 1.0E-6)

        flatCoordinates = fieldmodule.findFieldByName("flat coordinates").castFiniteElement()
        self.assertTrue(flatCoordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(flatCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [ 0.0, 0.0, 0.0 ], 1.0E-6)
        assertAlmostEqualList(self, maximums, [ 186.72988844629867, 77.4178187926561, 3.2000000000000006 ], 1.0E-6)

        textureCoordinates = fieldmodule.findFieldByName("texture coordinates").castFiniteElement()
        minimums, maximums = evaluateFieldNodesetRange(textureCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [ 0.0, 0.0, 0.0 ], 1.0E-6)
        assertAlmostEqualList(self, maximums, [ 0.9812471574796385, 1.0, 2.0 ], 1.0E-6)

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
        self.assertAlmostEqual(surfaceArea, 14612.416788520026, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 26826.069540921028, delta=1.0E-6)

    def test_mousecolon1(self):
        """
        Test creation of mouse colon scaffold.
        """
        options = MeshType_3d_colon1.getDefaultOptions("Mouse 2")
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = MeshType_3d_colon1.generateBaseMesh(region, options)
        self.assertEqual(6, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(400, mesh3d.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        flatCoordinates = fieldmodule.findFieldByName("flat coordinates").castFiniteElement()
        self.assertTrue(flatCoordinates.isValid())
        textureCoordinates = fieldmodule.findFieldByName("texture coordinates").castFiniteElement()
        self.assertTrue(textureCoordinates.isValid())

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            faceMeshGroup = createFaceMeshGroupExteriorOnFace(fieldmodule, Element.FACE_TYPE_XI3_1)
            flatSurfaceAreaField = fieldmodule.createFieldMeshIntegral(one, flatCoordinates, faceMeshGroup)
            flatSurfaceAreaField.setNumbersOfPoints(4)
            textureSurfaceAreaField = fieldmodule.createFieldMeshIntegral(one, textureCoordinates, faceMeshGroup)
            textureSurfaceAreaField.setNumbersOfPoints(4)
            textureVolumeField = fieldmodule.createFieldMeshIntegral(one, textureCoordinates, mesh3d)
            textureVolumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, flatSurfaceArea = flatSurfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(flatSurfaceArea, 652.2383326727861, delta=1.0E-6)
        result, textureSurfaceArea = textureSurfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(textureSurfaceArea, 1.0, delta=1.0E-6)
        result, textureVolume = textureVolumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(textureVolume, 1.0, delta=1.0E-6)

if __name__ == "__main__":
    unittest.main()
