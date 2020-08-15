'''
Class for generating a shield-shaped mesh, with a regular flat top but
a rounded bottom formed by having two points where 3 square elements
merge to form a triangle.
'''
from __future__ import division
import copy
import math
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.utils import vector, mirror
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, setEftScaleFactorIds
from scaffoldmaker.utils.interpolation import DerivativeScalingMode, sampleCubicHermiteCurves, smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition, calculate_surface_axes
from enum import Enum

class ShieldMode(Enum):
    SHIELD_MODE_FULL = 1        # full shield
    SHIELD_MODE_LOWER_HALF = 2  # half shield
    SHIELD_MODE_UPPER_HALF = 3  # half shield

class ShieldMesh:
    '''
    Shield mesh generator. Has one element through thickness.
    '''

    def __init__(self, elementsCountAcross, elementsCountUp, elementsCountRim, trackSurface : TrackSurface=None,
                 elementsCountAlong=1, shieldMode = ShieldMode.SHIELD_MODE_LOWER_HALF):
        '''
        Data structure for defining a shield-shaped mesh which is flat on the top and rounded around the bottom.
        It is represented as a regular box of elementsCountAcross x elementsCountUp
        but with strips of elements absent on the bottom left and right.
        Example showing elements for 4 across, 2 up and 0 rim, and how nodes/coordinates are stored in ShieldMesh:

        N___N___N___N___N           N___N___N___N___N
        |   |   |   |   |               |   |   |
        |   |   |   |   |               |   |   |
        \   T---N---T   /    -->        T---N---T
         \ /    |    \ /                |   |   |
          N     |     N                 |   |   |
           \----N----/                  N---N---N

        The removal of 2 corners of the box is achieved by triple points T where 3 square elements
        join in a triangle.
        Extra rim elements go around the sides and bottom curve. Rim elements add nodes on the sides above the
        triple points, and across the shorter bottom row.
        Extra elements up add regular rows of nodes/elements on top, and extra non-rim elements across
        add regular columns of nodes/elements up the centre.
        The mesh currently only supports 1 element through the thickness of the wall.
        :param elementsCountAcross: Number of elements across top of shield. Must be at least  4 + elementsCountRim.
        :param elementsCountUp: Number of elements up central axis of shield. Must be at least 2 + elementsCountRim.
        :param elementsCountRim: Number of elements around bottom rim (not top) outside of 'triple points'.
        :param trackSurface: Optional trackSurface to store or restrict points to.
        '''
        assert elementsCountRim >= 0
        assert elementsCountAlong >= 1
        assert elementsCountAcross >= (elementsCountRim + 4)
        assert elementsCountUp >= (elementsCountRim + 2)
        self.elementsCountAcross = elementsCountAcross
        self.elementsCountUp = elementsCountUp
        self.elementsCountRim = elementsCountRim
        self.elementsCountAlong = elementsCountAlong
        self.elementsCountUpRegular = elementsCountUp - 2 - elementsCountRim
        elementsCountAcrossNonRim = self.elementsCountAcross - 2*elementsCountRim
        self.elementsCountAroundFull = 2*self.elementsCountUpRegular + elementsCountAcrossNonRim
        self.trackSurface = trackSurface
        self._mode = shieldMode
        self.px  = [ [] for _ in range(elementsCountAlong+1) ]
        self.pd1 = [ [] for _ in range(elementsCountAlong+1) ]
        self.pd2 = [ [] for _ in range(elementsCountAlong+1) ]
        self.pd3 = [ [] for _ in range(elementsCountAlong+1) ]
        self.nodeId = [ [] for _ in range(elementsCountAlong+1) ]
        self._halvesCount = 2 if self._mode == ShieldMode.SHIELD_MODE_FULL else 1
        for n3 in range(elementsCountAlong+1):
            for n2 in range(self._halvesCount*elementsCountUp + 1):
                for p in [ self.px[n3], self.pd1[n3], self.pd2[n3], self.pd3[n3], self.nodeId[n3] ]:
                    p.append([ None ]*(elementsCountAcross + 1))
        if trackSurface:
            self.pProportions = [ [ None ]*(elementsCountAcross + 1) for n2 in range(elementsCountUp + 1) ]
        self.elementId = [ [ None ]*elementsCountAcross for n2 in range(self._halvesCount*elementsCountUp) ]


    def convertRimIndex(self, ix, rx=0):
        '''
        Convert point index around the lower rim to n1, n2 across and up box.
        :param ix: index around from 0 to self.elementsCountAroundFull
        :param rx: rim index from 0 (around outside) to self.elementsCountRim
        :return: n1, n2
        '''
        assert 0 <= ix <= self.elementsCountAroundFull
        assert 0 <= rx <= self.elementsCountRim
        if ix <= self.elementsCountUpRegular:
            return rx, self.elementsCountUp - ix
        mx = self.elementsCountAroundFull - ix
        if mx <= self.elementsCountUpRegular:
            return self.elementsCountAcross - rx, self.elementsCountUp - mx
        return self.elementsCountRim + ix - self.elementsCountUpRegular, rx


    def getTriplePoints(self, n3):
        '''
        Compute coordinates and derivatives of points where 3 square elements merge.
        :param n3: Index of through-wall coordinates to use.
        '''
        n1a = self.elementsCountRim
        n1b = n1a + 1
        n1c = n1a + 2
        m1a = self.elementsCountAcross - self.elementsCountRim
        m1b = m1a - 1
        m1c = m1a - 2
        n2a = self.elementsCountRim
        n2b = n2a + 1
        n2c = n2a + 2
        # left
        ltx = []
        tx, td1 = sampleCubicHermiteCurves(
            [ self.px[n3][n2a][n1c], self.px[n3][n2c][n1b] ], [ [ (-self.pd1[n3][n2a][n1c][c] - self.pd3[n3][n2a][n1c][c]) for c in range(3) ], self.pd1[n3][n2c][n1b] ], 2, arcLengthDerivatives = True)[0:2]
        ltx.append(tx[1])
        tx, td1 = sampleCubicHermiteCurves(
            [ self.px[n3][n2a][n1b], self.px[n3][n2c][n1c] ], [ [-self.pd3[n3][n2a][n1b][c] for c in range(3)], [ (self.pd1[n3][n2c][n1c][c] + self.pd3[n3][n2c][n1c][c]) for c in range(3) ] ], 2, arcLengthDerivatives = True)[0:2]
        ltx.append(tx[1])
        tx, td1 = sampleCubicHermiteCurves(
            [ self.px[n3][n2c][n1a], self.px[n3][n2b][n1c] ], [ [ (self.pd1[n3][n2c][n1a][c] - self.pd3[n3][n2c][n1a][c]) for c in range(3) ], self.pd3[n3][n2b][n1c] ], 2, arcLengthDerivatives = True)[0:2]
        ltx.append(tx[1])
        #x = [ (ltx[0][c] + ltx[1][c] + ltx[2][c])/3.0 for c in range(3) ]
        x = [ (ltx[0][c] + ltx[2][c])/2.0 for c in range(3) ]
        if self.trackSurface:
            p = self.trackSurface.findNearestPosition(x, startPosition=self.trackSurface.createPositionProportion(*(self.pProportions[n2b][n1c])))
            self.pProportions[n2b][n1b] = self.trackSurface.getProportion(p)
            x, sd1, sd2 = self.trackSurface.evaluateCoordinates(p, derivatives=True)
            d1, d2, d3 = calculate_surface_axes(sd1, sd2, vector.normalise(sd1))
            self.pd3[n3][n2b][n1b] = d3
        self.px [n3][n2b][n1b] = x
        self.pd3[n3][n2b][n1b] = [ (self.px[n3][n2b][n1c][c] - self.px[n3][n2b][n1b][c]) for c in range(3) ]
        self.pd1[n3][n2b][n1b] = [ (self.px[n3][n2c][n1b][c] - self.px[n3][n2b][n1b][c]) for c in range(3) ]
        if not self.trackSurface:
            self.pd2[n3][n2b][n1b] = vector.normalise(vector.crossproduct3(self.pd3[n3][n2b][n1b], self.pd1[n3][n2b][n1b]))
        # right
        rtx = []
        tx, td1 = sampleCubicHermiteCurves(
            [ self.px[n3][n2a][m1c], self.px[n3][n2c][m1b] ], [ [ (self.pd1[n3][n2a][m1c][c] - self.pd3[n3][n2a][m1c][c]) for c in range(3) ], self.pd1[n3][n2c][m1b] ], 2, arcLengthDerivatives = True)[0:2]
        rtx.append(tx[1])
        tx, td1 = sampleCubicHermiteCurves(
            [ self.px[n3][n2a][m1b], self.px[n3][n2c][m1c] ], [ [-self.pd3[n3][n2a][m1b][c] for c in range(3)], [ (-self.pd3[n3][n2c][m1c][c] + self.pd1[n3][n2c][m1c][c]) for c in range(3) ] ], 2, arcLengthDerivatives = True)[0:2]
        rtx.append(tx[1])
        tx, td1 = sampleCubicHermiteCurves(
            [ self.px[n3][n2c][m1a], self.px[n3][n2b][m1c] ], [ [ (-self.pd1[n3][n2c][m1a][c] - self.pd3[n3][n2c][m1a][c]) for c in range(3) ], [ -d for d in self.pd3[n3][n2b][m1c] ] ], 2, arcLengthDerivatives = True)[0:2]
        rtx.append(tx[1])
        #x = [ (rtx[0][c] + rtx[1][c] + rtx[2][c])/3.0 for c in range(3) ]
        x = [ (rtx[0][c] + rtx[2][c])/2.0 for c in range(3) ]
        if self.trackSurface:
            p = self.trackSurface.findNearestPosition(x, startPosition=self.trackSurface.createPositionProportion(*(self.pProportions[n2b][m1c])))
            self.pProportions[n2b][m1b] = self.trackSurface.getProportion(p)
            x, sd1, sd2 = self.trackSurface.evaluateCoordinates(p, derivatives=True)
            d1, d2, d3 = calculate_surface_axes(sd1, sd2, vector.normalise(sd1))
            self.pd3[n3][n2b][m1b] = d3
        self.px [n3][n2b][m1b] = x
        self.pd3[n3][n2b][m1b] = [ (self.px[n3][n2b][m1b][c] - self.px[n3][n2b][m1c][c]) for c in range(3) ]
        self.pd1[n3][n2b][m1b] = [ (self.px[n3][n2c][m1b][c] - self.px[n3][n2b][m1b][c]) for c in range(3) ]
        if not self.trackSurface:
            self.pd2[n3][n2b][m1b] = vector.normalise(vector.crossproduct3(self.pd3[n3][n2b][m1b], self.pd1[n3][n2b][m1b]))


    def smoothDerivativesToTriplePoints(self, n3, fixAllDirections=False):
        '''
        Smooth derivatives leading to triple points where 3 square elements merge.
        :param n3: Index of through-wall coordinates to use.
        '''
        n1a = self.elementsCountRim
        n1b = n1a + 1
        m1a = self.elementsCountAcross - self.elementsCountRim
        m1b = m1a - 1
        n2a = self.elementsCountRim
        n2b = n2a + 1
        n2c = n2a + 2
        # left
        tx = []
        td3 = []
        for n2 in range(0, n2c):
            tx .append(self.px [n3][n2][n1b])
            td3.append([-d for d in self.pd3[n3][n2][n1b]] if (n2 < n2b) else [ (self.pd1[n3][n2][n1b][c] + self.pd3[n3][n2][n1b][c]) for c in range(3) ])
        td3 = smoothCubicHermiteDerivativesLine(tx, td3, fixAllDirections=fixAllDirections, fixEndDerivative=True, magnitudeScalingMode = DerivativeScalingMode.HARMONIC_MEAN)
        for n2 in range(0, n2b):
            self.pd3[n3][n2][n1b] = [-d for d in td3[n2]]
        # right
        tx = []
        td3 = []
        for n2 in range(0, n2c):
            tx .append(self.px [n3][n2][m1b])
            td3.append([-d for d in self.pd3[n3][n2][m1b]] if (n2 < n2b) else [ (-self.pd3[n3][n2][m1b][c] + self.pd1[n3][n2][m1b][c]) for c in range(3) ])
        td3 = smoothCubicHermiteDerivativesLine(tx, td3, fixAllDirections=fixAllDirections, fixEndDerivative=True, magnitudeScalingMode = DerivativeScalingMode.HARMONIC_MEAN)
        for n2 in range(0, n2b):
            self.pd3[n3][n2][m1b] = [-d for d in td3[n2]]


    def smoothDerivativesAroundRim(self, n3, n3d=None, rx=0):
        '''
        Smooth derivatives around rim.
        :param n3: Index of through-wall coordinates to use.
        :param n3d: Which n3 index to copy initial derivatives from. If None, use n3.
        :param rx: rim index from 0 (around outside) to self.elementsCountRim
        '''
        assert 0 <= rx <= self.elementsCountRim
        if not n3d:
            n3d = n3
        tx = []
        td1 = []
        for ix in range(self.elementsCountAroundFull + 1):
            n1, n2 = self.convertRimIndex(ix, rx)
            tx.append(self.px[n3][n2][n1])
            if n2 > self.elementsCountRim:  # regular rows
                if n1 <= self.elementsCountRim:
                    td1.append([ -d for d in self.pd2[n3d][n2][n1] ])
                else:
                    td1.append(self.pd2[n3d][n2][n1])
            else:
                td1.append(self.pd1[n3d][n2][n1])
        td1 = smoothCubicHermiteDerivativesLine(tx, td1)
        for ix in range(self.elementsCountAroundFull + 1):
            n1, n2 = self.convertRimIndex(ix, rx)
            if n2 > self.elementsCountRim:  # regular rows
                if n1 <= self.elementsCountRim:
                    self.pd2[n3][n2][n1] = [ -d for d in td1[ix] ]
                else:
                    self.pd2[n3][n2][n1] = td1[ix]
            else:
                self.pd1[n3][n2][n1] = td1[ix]

    def generateNodesForOtherHalf(self, p):
        """
        Generates coordinates and derivativatives for the other half by mirroring them. It keeps the d1 dirrection.
        :param p: mirror plane ax+by+cz=d in form of [a,b,c,d]
        :return:
        """
        for n2 in range(self.elementsCountUp):
            for n3 in range(self.elementsCountAlong + 1):
                for n1 in range(self.elementsCountAcross + 1):
                    if self.px[n3][n2][n1]:
                        self.px[n3][2*self.elementsCountUp-n2][n1] = mirror.mirrorImageOfPoint(self.px[n3][n2][n1], p)
                        self.pd1[n3][2*self.elementsCountUp-n2][n1] = mirror.mirrorVector(self.pd1[n3][n2][n1], p,
                                                                      keepSign=False)
                        self.pd2[n3][2*self.elementsCountUp-n2][n1] = mirror.mirrorVector(self.pd2[n3][n2][n1], p)
                        self.pd3[n3][2*self.elementsCountUp-n2][n1] = mirror.mirrorVector(self.pd3[n3][n2][n1], p)

    def generateNodes(self, fieldmodule, coordinates, startNodeIdentifier,p):
        """
        Create shield nodes from coordinates.
        :param fieldmodule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
        :param coordinates: Coordinate field to define.
        :param startNodeIdentifier: First node identifier to use.
        :param p: mirror plane ax+by+cz=d in form of [a,b,c,d]
        :return: next nodeIdentifier.
         """
        nodeIdentifier = startNodeIdentifier
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        cache = fieldmodule.createFieldcache()

        #for n2 in range(self.elementsCountUp, -1, -1):
        #    s = ""
        #    for n1 in range(self.elementsCountAcross + 1):
        #        s += str(n1) if self.px[1][n2][n1] else " "
        #    print(n2, s, n2 - self.elementsCountUp - 1)

        if self._mode == ShieldMode.SHIELD_MODE_FULL:
            self.generateNodesForOtherHalf(p)

        for n2 in range(self._halvesCount*self.elementsCountUp + 1):
            for n3 in range(self.elementsCountAlong+1):
                for n1 in range(self.elementsCountAcross + 1):
                    if self.px[n3][n2][n1]:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        self.nodeId[n3][n2][n1] = nodeIdentifier
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, self.px [n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, self.pd1[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, self.pd2[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, self.pd3[n3][n2][n1])
                        nodeIdentifier += 1

        return nodeIdentifier


    def generateElements(self, fieldmodule, coordinates, startElementIdentifier, meshGroups=[]):
        """
        Create shield elements from nodes.
        :param fieldmodule: Zinc fieldmodule to create elements in.
        :param coordinates: Coordinate field to define.
        :param startElementIdentifier: First element identifier to use.
        :param meshGroups: Zinc mesh groups to add elements to.
        :return: next elementIdentifier.
         """
        elementIdentifier = startElementIdentifier
        useCrossDerivatives = False
        mesh = fieldmodule.findMeshByDimension(3)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftNoCrossDerivatives()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        isEven = (self.elementsCountAcross % 2) == 0
        e1a = self.elementsCountRim
        e1b = e1a + 1
        e1z = self.elementsCountAcross - 1 - self.elementsCountRim
        e1y = e1z - 1
        e2a = self.elementsCountRim
        e2b = self.elementsCountRim + 1
        e2c = self.elementsCountRim + 2
        e2d = 2*self.elementsCountUp-1
        for e3 in range(self.elementsCountAlong):
            for e2 in range(self._halvesCount*self.elementsCountUp):
                for e1 in range(self.elementsCountAcross):
                    eft1 = eft
                    scalefactors = None
                    nids = [ self.nodeId[e3][e2][e1], self.nodeId[e3][e2 + 1][e1], self.nodeId[e3+1][e2][e1], self.nodeId[e3+1][e2 + 1][e1],
                             self.nodeId[e3][e2][e1 + 1], self.nodeId[e3][e2 + 1][e1 + 1], self.nodeId[e3+1][e2][e1 + 1], self.nodeId[e3+1][e2 + 1][e1 + 1] ]
                    if (e2 < e2b) or (e2 == e2d):
                        if (e1 < e1b) or (e1 > e1y):
                            continue  # no element due to triple point closure
                        if (e2 == e2a) or (e2 == e2d):
                            # bottom and top row elements
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]
                            if e2 == e2a:
                                remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,[(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, [])])
                                if (e1 == e1b) or (e1 == e1y):
                                    # map bottom triple point element
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                    else:
                                        remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            elif e2 == e2d:
                                remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,[(Node.VALUE_LABEL_D_DS3, [])])
                                if (e1 == e1b) or (e1 == e1y):
                                    # map top triple point element
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS1,[(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1,[(Node.VALUE_LABEL_D_DS1, []),(Node.VALUE_LABEL_D_DS3, [])])

                    elif (e2 == e2b) or (e2 == e2d-e2b):
                        if (e1 <= e1a) or (e1 >= e1z):
                            # map top 2 triple point elements
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [ -1.0 ]
                            if e1 < e1a:
                                e2r = e1
                                nids[0] = self.nodeId[0][e2r    ][e1b]
                                nids[1] = self.nodeId[0][e2r + 1][e1b]
                                nids[4] = self.nodeId[1][e2r    ][e1b]
                                nids[5] = self.nodeId[1][e2r + 1][e1b]
                                remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                            elif e1 == e1a:
                                if e2 == e2b:
                                    nids[0] = self.nodeId[e3][e2a][e1b]
                                    nids[2] = self.nodeId[e3+1][e2a][e1b]
                                    tripleN = [5, 7]
                                    remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                elif e2 == e2d-e2b:
                                    nids[1] = self.nodeId[e3][e2d+1][e1b]
                                    nids[3] = self.nodeId[e3+1][e2d+1][e1b]
                                    tripleN = [6, 8]
                                    remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])

                            elif e1 == e1z:
                                if e2 == e2b:
                                    nids[4] = self.nodeId[e3][e2a][e1z]
                                    nids[6] = self.nodeId[e3+1][e2a][e1z]
                                    remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                elif e2 == e2d-e2b:
                                    nids[5] = self.nodeId[e3][e2d+1][e1z]
                                    nids[7] = self.nodeId[e3+1][e2d+1][e1z]
                                    remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                            elif e1 > e1z:
                                e2r = self.elementsCountAcross - e1
                                nids[0] = self.nodeId[0][e2r    ][e1z]
                                nids[1] = self.nodeId[0][e2r - 1][e1z]
                                nids[4] = self.nodeId[1][e2r    ][e1z]
                                nids[5] = self.nodeId[1][e2r - 1][e1z]
                                remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                    else:
                        if (e1 <= e1a):
                            # map left column elements
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [ -1.0 ]
                            remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS3, [1])])

                    if eft1 is not eft:
                        elementtemplate1.defineField(coordinates, -1, eft1)
                        element = mesh.createElement(elementIdentifier, elementtemplate1)
                    else:
                        element = mesh.createElement(elementIdentifier, elementtemplate)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    if scalefactors:
                        result3 = element.setScaleFactors(eft1, scalefactors)
                    else:
                        result3 = 7
                    #print('create element shield', elementIdentifier, result2, result3, nids)
                    self.elementId[e2][e1] = elementIdentifier
                    elementIdentifier += 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

        return elementIdentifier
