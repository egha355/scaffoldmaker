"""
Class for listing and accessing all mesh type scripts supported by scaffoldmaker.
"""

import json
from scaffoldmaker.meshtypes.meshtype_1d_bifurcationtree1 import MeshType_1d_bifurcationtree1
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.meshtypes.meshtype_2d_plate1 import MeshType_2d_plate1
from scaffoldmaker.meshtypes.meshtype_2d_platehole1 import MeshType_2d_platehole1
from scaffoldmaker.meshtypes.meshtype_2d_sphere1 import MeshType_2d_sphere1
from scaffoldmaker.meshtypes.meshtype_2d_tube1 import MeshType_2d_tube1
from scaffoldmaker.meshtypes.meshtype_2d_tubebifurcation1 import MeshType_2d_tubebifurcation1
from scaffoldmaker.meshtypes.meshtype_3d_bladder1 import MeshType_3d_bladder1
from scaffoldmaker.meshtypes.meshtype_3d_bladderurethra1 import MeshType_3d_bladderurethra1
from scaffoldmaker.meshtypes.meshtype_3d_box1 import MeshType_3d_box1
from scaffoldmaker.meshtypes.meshtype_3d_boxhole1 import MeshType_3d_boxhole1
from scaffoldmaker.meshtypes.meshtype_3d_cecum1 import MeshType_3d_cecum1
from scaffoldmaker.meshtypes.meshtype_3d_colon1 import MeshType_3d_colon1
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import MeshType_3d_colonsegment1
from scaffoldmaker.meshtypes.meshtype_3d_heart1 import MeshType_3d_heart1
from scaffoldmaker.meshtypes.meshtype_3d_heart2 import MeshType_3d_heart2
from scaffoldmaker.meshtypes.meshtype_3d_heartarterialroot1 import MeshType_3d_heartarterialroot1
from scaffoldmaker.meshtypes.meshtype_3d_heartarterialvalve1 import MeshType_3d_heartarterialvalve1
from scaffoldmaker.meshtypes.meshtype_3d_heartatria1 import MeshType_3d_heartatria1
from scaffoldmaker.meshtypes.meshtype_3d_heartatria2 import MeshType_3d_heartatria2
from scaffoldmaker.meshtypes.meshtype_3d_heartventricles1 import MeshType_3d_heartventricles1
from scaffoldmaker.meshtypes.meshtype_3d_heartventricles2 import MeshType_3d_heartventricles2
from scaffoldmaker.meshtypes.meshtype_3d_heartventricles3 import MeshType_3d_heartventricles3
from scaffoldmaker.meshtypes.meshtype_3d_heartventriclesbase1 import MeshType_3d_heartventriclesbase1
from scaffoldmaker.meshtypes.meshtype_3d_heartventriclesbase2 import MeshType_3d_heartventriclesbase2
from scaffoldmaker.meshtypes.meshtype_3d_lens1 import MeshType_3d_lens1
from scaffoldmaker.meshtypes.meshtype_3d_lung1 import MeshType_3d_lung1
from scaffoldmaker.meshtypes.meshtype_3d_ostium1 import MeshType_3d_ostium1
from scaffoldmaker.meshtypes.meshtype_3d_smallintestine1 import MeshType_3d_smallintestine1
from scaffoldmaker.meshtypes.meshtype_3d_solidsphere1 import MeshType_3d_solidsphere1
from scaffoldmaker.meshtypes.meshtype_3d_solidcylinder1 import MeshType_3d_solidcylinder1
from scaffoldmaker.meshtypes.meshtype_3d_sphereshell1 import MeshType_3d_sphereshell1
from scaffoldmaker.meshtypes.meshtype_3d_sphereshellseptum1 import MeshType_3d_sphereshellseptum1
from scaffoldmaker.meshtypes.meshtype_3d_stellate1 import MeshType_3d_stellate1
from scaffoldmaker.meshtypes.meshtype_3d_stomachhuman1 import MeshType_3d_stomachhuman1
from scaffoldmaker.meshtypes.meshtype_3d_tube1 import MeshType_3d_tube1
from scaffoldmaker.meshtypes.meshtype_3d_tubeseptum1 import MeshType_3d_tubeseptum1
from scaffoldmaker.meshtypes.meshtype_3d_wholebody1 import MeshType_3d_wholebody1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage


class Scaffolds(object):

    def __init__(self):
        self._allScaffoldTypes = [
            MeshType_1d_bifurcationtree1,
            MeshType_1d_path1,
            MeshType_2d_plate1,
            MeshType_2d_platehole1,
            MeshType_2d_sphere1,
            MeshType_2d_tube1,
            MeshType_2d_tubebifurcation1,
            MeshType_3d_bladder1,
            MeshType_3d_bladderurethra1,
            MeshType_3d_box1,
            MeshType_3d_boxhole1,
            MeshType_3d_cecum1,
            MeshType_3d_colon1,
            MeshType_3d_colonsegment1,
            MeshType_3d_heart1,
            MeshType_3d_heart2,
            MeshType_3d_heartarterialroot1,
            MeshType_3d_heartarterialvalve1,
            MeshType_3d_heartatria1,
            MeshType_3d_heartatria2,
            MeshType_3d_heartventricles1,
            MeshType_3d_heartventricles2,
            MeshType_3d_heartventricles3,
            MeshType_3d_heartventriclesbase1,
            MeshType_3d_heartventriclesbase2,
            MeshType_3d_lens1,
            MeshType_3d_lung1,
            MeshType_3d_ostium1,
            MeshType_3d_smallintestine1,
            MeshType_3d_solidsphere1,
            MeshType_3d_solidcylinder1,
            MeshType_3d_sphereshell1,
            MeshType_3d_sphereshellseptum1,
            MeshType_3d_stellate1,
            MeshType_3d_stomachhuman1,
            MeshType_3d_tube1,
            MeshType_3d_tubeseptum1,
            MeshType_3d_wholebody1
            ]

    def findScaffoldTypeByName(self, name):
        for scaffoldType in self._allScaffoldTypes:
            if scaffoldType.getName() == name:
                return scaffoldType
        return None

    def getDefaultMeshType(self):
        '''
        Deprecated: use getDefaultScaffoldType()
        '''
        return self.getDefaultScaffoldType()

    def getDefaultScaffoldType(self):
        return MeshType_3d_box1

    def getMeshTypes(self):
        '''
        Deprecated: use getScaffoldTypes()
        '''
        return self.getScaffoldTypes()

    def getScaffoldTypes(self):
        return self._allScaffoldTypes


class Scaffolds_JSONEncoder(json.JSONEncoder):
    '''
    Class encoding scaffold objects in JSON. Pass as cls argument to json.dumps.
    '''

    def default(self, obj):
        if isinstance(obj, ScaffoldPackage):
            dct = obj.toDict()
            dct['_ScaffoldPackage'] = True
            return dct
        elif isinstance(obj, bytes):
            return obj.decode("utf-8")
        else:
            super().default(obj)


def Scaffolds_decodeJSON(dct):
    '''
    Function for passing as object_hook argument to json.loads.
    Constructs scaffold objects from their JSON object encoding.
    '''
    if ('_ScaffoldPackage' in dct):
        scaffoldType = Scaffolds().findScaffoldTypeByName(dct['scaffoldTypeName'])
        #print('Scaffolds_decodeJSON scaffoldType',scaffoldType.getName(), dct)
        return ScaffoldPackage(scaffoldType, dct)
    return dct
