'''
Definitions of standard element field templates using bicubic Hermite x linear Lagrange basis.
'''
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, setEftScaleFactorIds
from opencmiss.zinc.element import Elementbasis, Elementfieldtemplate
from opencmiss.zinc.node import Node
from opencmiss.zinc.status import OK as ZINC_OK

class eftfactory_bicubichermitelinear:
    '''
    Factory class for creating element field templates for a 3-D mesh using bicubic Hermite x linear Lagrange basis.
    '''

    def __init__(self, mesh, useCrossDerivatives, linearAxis = 3,
            d_ds1 = Node.VALUE_LABEL_D_DS1, d_ds2 = Node.VALUE_LABEL_D_DS2):
        '''
        :param mesh:  Zinc mesh to create element field templates in.
        :param useCrossDerivatives: Set to True if you want cross derivative terms.
        :param linearAxis: 1, 2, or 3.
        :param d_ds1: Node derivative to use in Hermite axis 1: Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2.
        :param d_ds2: Node derivative to use in Hermite axis 2, > d_ds1: Node.VALUE_LABEL_D_DS2 or Node.VALUE_LABEL_D_DS3.
        '''
        assert mesh.getDimension() == 3, 'eftfactory_bicubichermitelinear: not a 3-D Zinc mesh'
        assert linearAxis in [ 1, 2, 3 ], 'eftfactory_bicubichermitelinear: linearAxis must be 1, 2 or 3'
        assert d_ds1 in [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2 ], 'eftfactory_bicubichermitelinear: invalid d_ds1'
        assert d_ds2 in [ Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ] and (d_ds2 > d_ds1), 'eftfactory_bicubichermitelinear: invalid d_ds2'
        self._mesh = mesh
        self._useCrossDerivatives = useCrossDerivatives
        self._linearAxis = linearAxis
        self._d_ds1 = d_ds1
        self._d_ds2 = d_ds2
        self._d2_ds1ds2 = Node.VALUE_LABEL_D2_DS2DS3 if (d_ds1 == Node.VALUE_LABEL_D_DS2) \
                     else Node.VALUE_LABEL_D2_DS1DS3 if (d_ds2 == Node.VALUE_LABEL_D_DS3) \
                     else Node.VALUE_LABEL_D2_DS1DS2
        self._fieldmodule = mesh.getFieldmodule()
        self._basis = self._fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        self._basis.setFunctionType(linearAxis, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)

    def _remapDefaultNodeDerivatives(self, eft):
        '''
        Remap the Hermite node derivatives to those chosen in __init__.
        Use only on first create.
        :param eft: The element field template to remap.
        '''
        # must do d_ds2 first!
        if self._d_ds2 != Node.VALUE_LABEL_D_DS2:
            remapEftNodeValueLabel(eft, range(1, 9), Node.VALUE_LABEL_D_DS2, [ (self._d_ds2, []) ])
        if self._d_ds1 != Node.VALUE_LABEL_D_DS1:
            remapEftNodeValueLabel(eft, range(1, 9), Node.VALUE_LABEL_D_DS1, [ (self._d_ds1, []) ])
        if self._d2_ds1ds2 != Node.VALUE_LABEL_D2_DS1DS2:
            remapEftNodeValueLabel(eft, range(1, 9), Node.VALUE_LABEL_D2_DS1DS2, [ (self._d2_ds1ds2, []) ])

    def createEftBasic(self):
        '''
        Create the basic biicubic Hermite x linear Lagrange element template with 1:1 mappings to
        node derivatives ds1 & ds2, with or without cross derivatives accordinate as initialised.
        :return: Element field template
        '''
        if not self._useCrossDerivatives:
            return self.createEftNoCrossDerivatives()
        eft = self._mesh.createElementfieldtemplate(self._basis)
        self._remapDefaultNodeDerivatives(eft)
        assert eft.validate(), 'eftfactory_bicubichermitelinear.createEftBasic:  Failed to validate eft'
        return eft

    def createEftNoCrossDerivatives(self):
        '''
        Create a basic tricubic hermite element template with 1:1 mappings to
        node derivatives ds1 & ds2, without cross derivatives.
        :return: Element field template
        '''
        eft = self._mesh.createElementfieldtemplate(self._basis)
        for n in range(8):
            eft.setFunctionNumberOfTerms(n*4 + 4, 0)
        self._remapDefaultNodeDerivatives(eft)
        assert eft.validate(), 'eftfactory_bicubichermitelinear.createEftNoCrossDerivatives:  Failed to validate eft'
        return eft

    def createEftShellPoleBottom(self, nodeScaleFactorOffset0, nodeScaleFactorOffset1):
        '''
        Create a bicubic hermite linear element field for closing bottom pole of a shell.
        Element is collapsed in xi1 on xi2 = 0.
        Each collapsed node has 3 scale factors giving the cos, sin coefficients
        of the radial line from global derivatives, plus the arc subtended by
        the element in radians, so the pole can be rounded.
        Need to create a new template for each sector around pole giving common
        nodeScaleFactorOffset values on common faces. Suggestion is to start at 0 and
        add 100 for each radial line around pole.
        :param nodeScaleFactorOffset0: offset of node scale factors at pole on xi1=0
        :param nodeScaleFactorOffset1: offset of node scale factors at pole on xi1=1
        :return: Element field template
        '''
        # start with full bicubic hermite linear to remap D2_DS1DS2 at pole
        eft = self._mesh.createElementfieldtemplate(self._basis)
        if not self._useCrossDerivatives:
            for n in [ 2, 3, 6, 7 ]:
                eft.setFunctionNumberOfTerms(n*4 + 4, 0)

        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        setEftScaleFactorIds(eft, [1], [
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3,
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3 ])
        # remap parameters before collapsing nodes
        remapEftNodeValueLabel(eft, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [])
        for layer in range(2):
            so = layer*6 + 1
            ln = layer*4 + 1
            # 2 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [so + 1]), (Node.VALUE_LABEL_D_DS2, [so + 2]) ])
            # 2 terms for cross derivative 1 2 to correct circular pole: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [so + 2, so + 3]), (Node.VALUE_LABEL_D_DS2, [1, so + 1, so + 3]) ])

            ln = layer*4 + 2
            # 2 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [so + 4]), (Node.VALUE_LABEL_D_DS2, [so + 5]) ])
            # 2 terms for cross derivative 1 2 to correct circular pole: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [so + 5, so + 6]), (Node.VALUE_LABEL_D_DS2, [1, so + 4, so + 6]) ])

        ln_map = [ 1, 1, 2, 3, 4, 4, 5, 6 ]
        remapEftLocalNodes(eft, 6, ln_map)

        assert eft.validate(), 'eftfactory_tricubichermite.createEftShellPoleBottom:  Failed to validate eft'
        return eft

    def createEftShellPoleTop(self, nodeScaleFactorOffset0, nodeScaleFactorOffset1):
        '''
        Create a bicubic hermite linear element field for closing top pole of a shell.
        Element is collapsed in xi1 on xi2 = 1.
        Each collapsed node has 3 scale factors giving the cos, sin coefficients
        of the radial line from global derivatives, plus the arc subtended by
        the element in radians, so the pole can be rounded.
        Need to create a new template for each sector around pole giving common
        nodeScaleFactorOffset values on common faces. Suggestion is to start at 0 and
        add 100 for each radial line around pole.
        :param nodeScaleFactorOffset0: offset of node scale factors at pole on xi1=0
        :param nodeScaleFactorOffset1: offset of node scale factors at pole on xi1=1
        :return: Element field template
        '''
        # start with full bicubic hermite linear to remap D2_DS1DS2 at pole
        eft = self._mesh.createElementfieldtemplate(self._basis)
        if not self._useCrossDerivatives:
            for n in [ 0, 1, 4, 5 ]:
                eft.setFunctionNumberOfTerms(n*4 + 4, 0)

        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        setEftScaleFactorIds(eft, [1], [
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3,
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3 ])
        # remap parameters before collapsing nodes
        remapEftNodeValueLabel(eft, [ 3, 4, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
        for layer in range(2):
            so = layer*6 + 1
            ln = layer*4 + 3
            # 2 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [so + 1]), (Node.VALUE_LABEL_D_DS2, [so + 2]) ])
            # 2 terms for cross derivative 1 2 to correct circular pole: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [1, so + 2, so + 3]), (Node.VALUE_LABEL_D_DS2, [so + 1, so + 3]) ])

            ln = layer*4 + 4
            # 2 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, so + 4), (Node.VALUE_LABEL_D_DS2, so + 5) ])
            # 2 terms for cross derivative 1 2 to correct circular pole: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [1, so + 5, so + 6]), (Node.VALUE_LABEL_D_DS2, [so + 4, so + 6]) ])

        ln_map = [ 1, 2, 3, 3, 4, 5, 6, 6 ]
        remapEftLocalNodes(eft, 6, ln_map)

        assert eft.validate(), 'eftfactory_tricubichermite.createEftShellPoleTop:  Failed to validate eft'
        return eft
        
    def createEftSplitXi1RightStraight(self):
        '''
        Create an element field template suitable for the inner elements of the
        join between left and right chambers, with xi1 bifurcating to right.
        Straight through version.
        Only works with linearAxis 2.
        :return: Element field template
        '''
        assert linearAxis == 2, 'eftfactory_bicubichermitelinear.createEftSplitXi1RightStraight:  Not linearAxis 2'
        eft = self.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eft, [1], [])
        remapEftNodeValueLabel(eft, [ 5, 7 ], self._d_ds1, [ (self._d_ds1, []), (self._d_ds2, [1]) ])
        remapEftNodeValueLabel(eft, [ 6, 8 ], self._d_ds1, [ (self._d_ds1, []), (self._d_ds2, []) ])
        assert eft.validate(), 'eftfactory_bicubichermitelinear.createEftSplitXi1RightStraight:  Failed to validate eft'
        return eft

    def createEftSplitXi1RightOut(self):
        '''
        Create an element field template suitable for the outer elements of the
        join between left and right chambers, with xi1 bifurcating to right.
        Right out version i.e. xi1 heading to right. h-shape.
        Only works with linearAxis 2.
        :return: Element field template
        '''
        assert linearAxis == 2, 'eftfactory_bicubichermitelinear.createEftSplitXi1RightOut:  Not linearAxis 2'
        eft = self.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eft, [1], [])
        remapEftNodeValueLabel(eft, [ 1, 3 ], self._d_ds1, [ (self._d_ds1, [1]) ])
        remapEftNodeValueLabel(eft, [ 1, 3 ], self._d_ds2, [ (self._d_ds1, [1]), (self._d_ds2, [1]) ])
        remapEftNodeValueLabel(eft, [ 5, 7 ], self._d_ds2, [ (self._d_ds1, [1]), (self._d_ds2, []) ])
        assert eft.validate(), 'eftfactory_bicubichermitelinear.createEftSplitXi1RightOut:  Failed to validate eft'
        return eft

    def createEftOpenTube(self):
        '''
        Create a basic bicubic hermite linear element template for elements
        along boundary where a tube is opened on xi1 = 1 for a flat preparation.
        Could eventually have 6 variants. Retain node numbering with two versions
        for boundary nodes.
        :return: Element field template
        '''
        eft = self.createEftBasic()
        for n in [ 1, 3, 5, 7 ]:
            ln = n + 1
            eft.setTermNodeParameter(n*4 + 1, 1, ln, Node.VALUE_LABEL_VALUE, 2)
            eft.setTermNodeParameter(n*4 + 2, 1, ln, Node.VALUE_LABEL_D_DS1, 2)
            eft.setTermNodeParameter(n*4 + 3, 1, ln, Node.VALUE_LABEL_D_DS2, 2)
            if self._useCrossDerivatives:
                eft.setTermNodeParameter(n*4 + 4, 1, ln, Node.VALUE_LABEL_D2_DS1DS2, 2)

        assert eft.validate(), 'eftfactory_bicubichermitelinear.createEftOpenTube:  Failed to validate eft'
        return eft

    def createEftWedgeXi1One(self):
        '''
        Create a basic bicubic hermite linear element template for elements
        along boundary of tenia coli where nodes on xi1 = 1 are collapsed.
        :return: Element field template
        '''
        eft = self.createEftBasic()
        ln_map = [ 1, 2, 3, 4, 5, 2, 6, 4 ]
        remapEftLocalNodes(eft, 6, ln_map)
        assert eft.validate(), 'eftfactory_tricubichermite.createEftWedgeXi1One:  Failed to validate eft'
        return eft

    def createEftWedgeXi1Zero(self):
        '''
        Create a basic bicubic hermite linear element template for elements
        along boundary of tenia coli where nodes on xi1 = 0 are collapsed.
        :return: Element field template
        '''
        eft = self.createEftBasic()
        ln_map = [ 1, 2, 3, 4, 1, 5, 3, 6 ]
        remapEftLocalNodes(eft, 6, ln_map)
        assert eft.validate(), 'eftfactory_tricubichermite.createEftWedgeXi1Zero:  Failed to validate eft'
        return eft

    def createEftWedgeCollapseXi1Quadrant(self, collapseNodes):
        '''
        Create a bicubic hermite linear element field for a wedge element collapsed in xi1.
        :param collapseNodes: As the element can be collapsed in xi1 at either ends of xi2 or xi3, collapseNodes
        are the local indices of nodes whose d2 (for elements collapse at either ends of xi2) or
        d3 (for elements collapse at either ends of xi3) are remapped with d1 before collapsing the nodes.
        :return: Element field template
        '''
        eft = self.createEftBasic()
        
        valid = True

        if collapseNodes in [[1, 3], [2, 4]]: # xi3 = 0
            nodes = [1, 2, 3, 4]
            remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
            ln_map = [1, 1, 2, 2, 3, 4, 5, 6]
        elif collapseNodes in [[5, 7], [6, 8]]: # xi3 = 1
            nodes = [5, 6, 7, 8]
            remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
            ln_map = [1, 2, 3, 4, 5, 5, 6, 6]
        elif collapseNodes in [[1, 5], [2, 6]]:
            nodes = [1, 2, 5, 6]
            # remap parameters on xi2 = 0 before collapsing nodes
            if collapseNodes == [1, 5]:
                setEftScaleFactorIds(eft, [1], [])
                remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
            elif collapseNodes == [2, 6]:
                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
            else:
                valid = False
            ln_map = [1, 1, 2, 3, 4, 4, 5, 6]
        elif collapseNodes in [[3, 7], [4, 8]]:
            nodes = [3, 4, 7, 8]
            # remap parameters on xi2 = 1 before collapsing nodes
            if collapseNodes == [3, 7]:
                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
            elif collapseNodes == [4, 8]:
                setEftScaleFactorIds(eft, [1], [])
                remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
            else:
                valid = False
            ln_map = [1, 2, 3, 3, 4, 5, 6, 6]
        else:
            valid = False

        if not valid:
            assert False, "createEftWedgeCollapseXi1Quadrant.  Not implemented for collapse nodes " + str(collapseNodes)

        # zero cross derivative parameters
        remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D2_DS1DS2, [])

        remapEftLocalNodes(eft, 6, ln_map)
        if not eft.validate():
            print('eftfactory_bicubichermitelinear.createEftWedgeCollapseXi1Quadrant:  Failed to validate eft for collapseNodes', collapseNodes)
        return eft

    def createEftWedgeCollapseXi2Quadrant(self, collapseNodes):
        '''
        Create a bicubic hermite linear element field for a wedge element collapsed in xi2.
        :param collapseNodes: As the element can be collapsed in xi2 at either ends of xi1 or xi3, collapseNodes
        are the local indices of nodes whose d1 (for elements collapse at either ends of xi1) or
        d3 (for elements collapse at either ends of xi3) are remapped with d2 before collapsing the nodes.
        :return: Element field template
        '''
        eft = self.createEftBasic()

        valid = True
        if collapseNodes in [[1, 2], [3, 4]]: # xi3 = 0
            nodes = [1, 2, 3, 4]
            remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS2, [])
            ln_map = [1, 2, 1, 2, 3, 4, 5, 6]
        elif collapseNodes in [[5, 6], [7, 8]]: # xi3 = 1
            nodes = [5, 6, 7, 8]
            remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS2, [])
            ln_map = [1, 2, 3, 4, 5, 6, 5, 6]

        elif collapseNodes in [[3, 7]]:
            nodes = [1, 3, 5, 7]
            # remap parameters on xi1 = 0 before collapsing nodes
            if collapseNodes == [3, 7]:
                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS2, [])
                remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
            else:
                valid = False
            ln_map = [1, 2, 1, 3, 4, 5, 4, 6]

        elif collapseNodes in [[2, 6], [4, 8]]:
            nodes = [2, 4, 6, 8]
            # remap parameters on xi1 = 1 before collapsing nodes
            if collapseNodes == [2, 6]:
                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS2, [])
                remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
            elif collapseNodes == [4, 8]:
                setEftScaleFactorIds(eft, [1], [])
                remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS2, [])
            else:
                valid = False
            ln_map = [1, 2, 3, 2, 4, 5, 6, 5]

        else:
            valid = False

        if not valid:
            assert False, "createEftWedgeCollapseXi2Quadrant.  Not implemented for collapse nodes " + str(collapseNodes)

        # zero cross derivative parameters
        remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D2_DS1DS2, [])

        remapEftLocalNodes(eft, 6, ln_map)
        if not eft.validate():
            print('eftfactory_bicubichermitelinear.createEftWedgeCollapseXi2Quadrant:  Failed to validate eft for collapseNodes', collapseNodes)
        return eft

    def createEftWedgeXi1ZeroOpenTube(self):
        '''
        Create a basic bicubic hermite linear element template for elements
        along boundary of tenia coli where nodes on xi1 = 0 are collapsed
        where a tube is opened on xi1 = 1 for a flat preparation.
        :return: Element field template
        '''
        eft = self.createEftBasic()
        for n in [ 1, 3, 5, 7 ]:
            ln = n + 1
            eft.setTermNodeParameter(n*4 + 1, 1, ln, Node.VALUE_LABEL_VALUE, 2)
            eft.setTermNodeParameter(n*4 + 2, 1, ln, Node.VALUE_LABEL_D_DS1, 2)
            eft.setTermNodeParameter(n*4 + 3, 1, ln, Node.VALUE_LABEL_D_DS2, 2)
            if self._useCrossDerivatives:
                eft.setTermNodeParameter(n*4 + 4, 1, ln, Node.VALUE_LABEL_D2_DS1DS2, 2)
        ln_map = [ 1, 2, 3, 4, 1, 5, 3, 6 ]
        remapEftLocalNodes(eft, 6, ln_map)
        assert eft.validate(), 'eftfactory_tricubichermite.createEftWedgeXi1ZeroOpenTube:  Failed to validate eft'
        return eft

    def createEftTetrahedronXi1One(self, nodeScaleFactorOffset0, nodeScaleFactorOffset1):
        '''
        Create a bicubic hermite linear element field for a solid tetrahedron for the apex of cecum,
        with xi1 and xi3 collapsed on xi2 = 0, and xi3 collapsed on xi1 = 1 and xi2 = 1.
        Each collapsed node on xi2 = 0 has 3 scale factors giving the cos, sin coefficients of
        the radial line from global derivatives, plus the arc subtended by the element in radians,
        so the circumferential direction is rounded.
        Need to create a new template for each sector around axis giving common
        nodeScaleFactorOffset values on common faces. Suggestion is to start at 0 and
        add 10000 for each radial line around axis.
        :param nodeScaleFactorOffset0: offset of node scale factors at axis on xi1=0
        :param nodeScaleFactorOffset1: offset of node scale factors at axis on xi1=1
        :return: Element field template
        '''
        # start with full bicubic hermite linear
        eft = self._mesh.createElementfieldtemplate(self._basis)

        for n in [ 2, 3, 6, 7 ]:
            eft.setFunctionNumberOfTerms(n * 4 + 4, 0)

        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        setEftScaleFactorIds(eft, [1], [
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3 ] )

        # remap parameters on xi2 = 0 before collapsing nodes
        remapEftNodeValueLabel(eft, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [])

        for layer in range(2):
            soAround = 1
            ln = layer * 4 + 1
            # 2 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ln], Node.VALUE_LABEL_D_DS2,
                                   [(Node.VALUE_LABEL_D_DS1, [soAround + 1]),
                                   (Node.VALUE_LABEL_D_DS2, [soAround + 2])])
            # 2 terms for cross derivative 1 2 to correct circular apex: cos(theta).phi, -sin(theta).phi
            remapEftNodeValueLabel(eft, [ln], Node.VALUE_LABEL_D2_DS1DS2,
                                   [(Node.VALUE_LABEL_D_DS1, [soAround + 2, soAround + 3]),
                                   (Node.VALUE_LABEL_D_DS2, [1, soAround + 1, soAround + 3])])

            ln = layer * 4 + 2
            # 2 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ln], Node.VALUE_LABEL_D_DS2,
                                   [(Node.VALUE_LABEL_D_DS1, [soAround + 4]),
                                   (Node.VALUE_LABEL_D_DS2, [soAround + 5])])
            # 2 terms for cross derivative 1 2 to correct circular apex: cos(theta).phi, -sin(theta).phi
            remapEftNodeValueLabel(eft, [ln], Node.VALUE_LABEL_D2_DS1DS2,
                                   [(Node.VALUE_LABEL_D_DS1, [soAround + 5, soAround + 6]),
                                    (Node.VALUE_LABEL_D_DS2, [1, soAround + 4, soAround + 6])])

        ln_map = [ 1, 1, 2, 3, 1, 1, 4, 3]
        remapEftLocalNodes(eft, 4, ln_map)

        assert eft.validate(), 'eftfactory_bicubichermitelinear.createEftTetrahedronXi1One:  Failed to validate eft'
        return eft

    def createEftTetrahedronXi1Zero(self, nodeScaleFactorOffset0, nodeScaleFactorOffset1):
        '''
        Create a bicubic hermite linear element field for a solid tetrahedron for the apex of cecum,
        with xi1 and xi3 collapsed on xi2 = 0, and xi3 collapsed on xi1 = 0, xi2 = 1.
        Each collapsed node on xi2 = 0 has 3 scale factors giving the cos, sin coefficients of
        the radial line from global derivatives, plus the arc subtended by the element in radians,
        so the circumferential direction is rounded.
        Need to create a new template for each sector around axis giving common
        nodeScaleFactorOffset values on common faces. Suggestion is to start at 0 and
        add 10000 for each radial line around axis.
        :param nodeScaleFactorOffset0: offset of node scale factors at axis on xi1=0
        :param nodeScaleFactorOffset1: offset of node scale factors at axis on xi1=1
        :return: Element field template
        '''
        # start with full bicubic hermite linear
        eft = self._mesh.createElementfieldtemplate(self._basis)
        for n in [ 2, 3, 6, 7 ]:
            eft.setFunctionNumberOfTerms(n * 4 + 4, 0)

        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        setEftScaleFactorIds(eft, [1], [
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3 ])

        # remap parameters on xi2 = 0 before collapsing nodes
        remapEftNodeValueLabel(eft, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [])

        for layer in range(2):
            soAround = 1
            ln = layer * 4 + 1
            # 2 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ln], Node.VALUE_LABEL_D_DS2,
                                   [(Node.VALUE_LABEL_D_DS1, [soAround + 1]),
                                   (Node.VALUE_LABEL_D_DS2, [soAround + 2])])
            # 2 terms for cross derivative 1 2 to correct circular apex: cos(theta).phi, -sin(theta).phi
            remapEftNodeValueLabel(eft, [ln], Node.VALUE_LABEL_D2_DS1DS2,
                                   [(Node.VALUE_LABEL_D_DS1, [soAround + 2, soAround + 3]),
                                   (Node.VALUE_LABEL_D_DS2, [1, soAround + 1, soAround + 3])])

            ln = layer * 4 + 2
            # 2 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ln], Node.VALUE_LABEL_D_DS2,
                                   [(Node.VALUE_LABEL_D_DS1, [soAround + 4]),
                                   (Node.VALUE_LABEL_D_DS2, [soAround + 5])])
            # 2 terms for cross derivative 1 2 to correct circular apex: cos(theta).phi, -sin(theta).phi
            remapEftNodeValueLabel(eft, [ln], Node.VALUE_LABEL_D2_DS1DS2,
                                   [(Node.VALUE_LABEL_D_DS1, [soAround + 5, soAround + 6]),
                                    (Node.VALUE_LABEL_D_DS2, [1, soAround + 4, soAround + 6])])

        ln_map = [ 1, 1, 2, 3, 1, 1, 2, 4]
        remapEftLocalNodes(eft, 4, ln_map)

        assert eft.validate(), 'eftfactory_bicubichermitelinear.createEftTetrahedronXi1Zero:  Failed to validate eft'
        return eft

    def createEftPyramidBottomSimple(self, nodeScaleFactorOffset0, nodeScaleFactorOffset1):
        '''
        Create a bicubic hermite linear element field for a solid pyramid for elements within
        a tenia coli joining to the cecal apex, with xi1 and xi3 collapsed on xi2 = 0.
        Each collapsed node has 3 scale factors giving the cos, sin coefficients of the
        radial line from global derivatives, plus the arc subtended by the element in radians,
        so the circumferential direction is rounded. Need to create a new template for each
        sector around axis giving common nodeScaleFactorOffset values on common faces.
        Suggestion is to start at 0 and add 10000 for each radial line around axis.
        :param nodeScaleFactorOffset0: offset of node scale factors at axis on xi1=0
        :param nodeScaleFactorOffset1: offset of node scale factors at axis on xi1=1
        :return: Element field template
        '''
        # start with full bicubic hermite linear
        eft = self._mesh.createElementfieldtemplate(self._basis)
        for n in [ 2, 3, 6, 7 ]:
            eft.setFunctionNumberOfTerms(n * 4 + 4, 0)

        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        setEftScaleFactorIds(eft, [1], [
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3])

        # remap parameters on xi2 = 0 before collapsing nodes
        remapEftNodeValueLabel(eft, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [])

        for layer in range(2):
            soAround = 1
            ln = layer * 4 + 1
            # 2 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2,
                                   [(Node.VALUE_LABEL_D_DS1, [soAround + 1]), (Node.VALUE_LABEL_D_DS2, [soAround + 2])])
            # 2 terms for cross derivative 1 2 to correct circular apex: cos(theta).phi, -sin(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2,
                                   [(Node.VALUE_LABEL_D_DS1, [soAround + 2, soAround + 3]),
                                    (Node.VALUE_LABEL_D_DS2, [1, soAround + 1, soAround + 3])])

            ln = layer * 4 + 2
            # 2 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2,
                                   [(Node.VALUE_LABEL_D_DS1, [soAround + 4]), (Node.VALUE_LABEL_D_DS2, [soAround + 5])])
            # 2 terms for cross derivative 1 2 to correct circular apex: cos(theta).phi, -sin(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2,
                                   [(Node.VALUE_LABEL_D_DS1, [soAround + 5, soAround + 6]),
                                    (Node.VALUE_LABEL_D_DS2, [1, soAround + 4, soAround + 6])])

        ln_map = [ 1, 1, 2, 3, 1, 1, 4, 5 ]
        remapEftLocalNodes(eft, 5, ln_map)

        assert eft.validate(), 'eftfactory_bicubichermitelinear.createEftPyramidBottomSimple:  Failed to validate eft'
        return eft
