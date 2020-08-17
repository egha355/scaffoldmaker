"""
Utility functions for generating a generalised 3-D solid cylinder (extruded ellipse/circle). It can be used to genrate a
solid truncated cone. It also can be used for transition from a 2D base to another base (e.g., ellipse to a circle).
"""

from enum import Enum
from scaffoldmaker.utils import vector, geometry, mirror
import math
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier, getMaximumElementIdentifier
from scaffoldmaker.utils.shieldmesh2 import ShieldMesh, ShieldMode
from scaffoldmaker.utils.interpolation import sampleCubicHermiteCurves, interpolateSampleCubicHermite,\
    smoothCubicHermiteDerivativesLine

class CylinderMode(Enum):
    CYLINDER_MODE_FULL = 1  # full cylinder is createted
    CYLINDER_MODE_LOWER_HALF = 2  # lower half cylinder

class CylinderType(Enum):
    CYLIDNER_REGULAR = 1        # all the bases along the axis of the cylinder are the same.
    CYLIDNER_TRUNCATED_CONE = 2 # different ellipses along the cylinder axis

class CylinderMesh:
    '''
    Cylinder mesh generator. Extrudes an ellipse/circle.
    '''

    def __init__(self,fieldModule, coordinates, baseCentre, alongAxis, majorAxis, minorRadius,
                             elementsCountAcross, elementsCountUp, elementsCountAlong,
                             cylinderMode = CylinderMode.CYLINDER_MODE_FULL,
                             cylinderType = CylinderType.CYLIDNER_REGULAR,
                             rate = 0.0, useCrossDerivatives = False):
        '''
        :param fieldmodule: Zinc fieldmodule to create elements in.
        :param coordinates: Coordinate field to define.
        :param baseCentre: The centre of the base of the cylinder (ellipse/circle).
        :param alongAxis: The cylinder axis that the base is extruded along.
        :param majorAxis: The major axis of the base. Should be perpendicular to alongAxis
        :param minorRadius: The minor radius of the base.
        :param elementsCountAcross: Number of elements across the base minor radius.
        :param elementsCountUp: Number of elements up central axis of shield. Must be at least 2 + elementsCountRim.
        :param elementsCountAlong: Number of elements along the cylinder axis.
        :param rate: Cone base radius reduction rate along cone axis.
        :param cylinderMode: A value from enum CylinderMode specifying.
        '''
        self._baseCentre = baseCentre
        self._alongAxis = alongAxis
        self._majorAxis = majorAxis
        self._minorRadius = minorRadius
        self._shield = None
        self._elementsCountAcross = elementsCountAcross
        self._elementsCountUp = elementsCountUp
        self._elementsCountAlong = elementsCountAlong
        self._startNodeIdentifier = 1
        self._startElementIdentifier = 1
        self._endNodeIdentifier = 1
        self._endElementIdentifier = 1
        self._cylinderMode = cylinderMode
        self._cylinderType = cylinderType
        self._useCrossDerivatives = useCrossDerivatives
        self._height = vector.magnitude(alongAxis)
        self._majorRadius = vector.magnitude(majorAxis)
        self._radiusReductionRate = rate
        # generate the mesh
        self.createCylinderMesh3d(fieldModule, coordinates)

    @staticmethod
    def createCylinderBaseMesh2D(centre, majorAxis, minorAxis, elementsCountAround, height):
        '''
        Generate a set of points and derivatives for an ellipse
        starting at pole majorAxis from centre.
        :param centre: Centre of full ellipse.
        :param majorAxis: Vector in direction of starting major radius, magnitude is ellipse major radius.
        :param minorAxis: Vector normal to major axis, magnitude is ellipse minor axis length.
        :param height: Height of arc of ellipsoid from starting point along majorAxis.
        :return: Lists nx, nd1. Ordered fastest around, starting at major radius.
        '''
        nx = []
        nd1 = []
        magMajorAxis = vector.magnitude(majorAxis)
        magMinorAxis = vector.magnitude(minorAxis)
        unitMajorAxis = vector.normalise(majorAxis)
        unitMinorAxis = vector.normalise(minorAxis)
        useHeight = min(max(0.0, height), 2.0 * magMajorAxis)
        totalRadians = geometry.getEllipseRadiansToX(magMajorAxis, 0.0, magMajorAxis - useHeight,
                                                     initialTheta=0.5 * math.pi * useHeight / magMajorAxis)
        radians = 0.0
        arcLengthUp = geometry.getEllipseArcLength(magMajorAxis, magMinorAxis, radians, totalRadians)
        elementsCountUp = elementsCountAround // 2
        elementArcLength = arcLengthUp / elementsCountUp
        radians = geometry.updateEllipseAngleByArcLength(magMajorAxis, magMinorAxis, radians, -arcLengthUp)
        for n1 in range(2*elementsCountUp + 1):
            cosRadians = math.cos(radians)
            sinRadians = math.sin(radians)
            nx.append(
                [(centre[c] + cosRadians * majorAxis[c] + sinRadians * minorAxis[c]) for c in range(3)])

            ndab = vector.setMagnitude([-sinRadians * magMajorAxis, cosRadians * magMinorAxis], elementArcLength)
            nd1.append(
                [(ndab[0] * unitMajorAxis[c] + ndab[1] * unitMinorAxis[c]) for c in range(3)])
            radians = geometry.updateEllipseAngleByArcLength(magMajorAxis, magMinorAxis, radians, elementArcLength)
        return nx, nd1

    def generateBasesMesh(self,majorRadius,elementsCountAround,arcLengthAlong,minorAxis,rate=0.08):
        '''
        generate bases of the truncated cone along the cone axis.
        :param majorRadius: major radius of the cone ellipse base.
        :param elementsCountAround: major radius of the cone ellipse base.
        :return:
        '''
        nx, nd1 = self.createCylinderBaseMesh2D(
            self._baseCentre, self._majorAxis, minorAxis, elementsCountAround, majorRadius)
        majorRadius1 = majorRadius
        tnx, tnd1, tnd2, tnd3 = [], [], [], []
        for n3 in range(self._elementsCountUp + 1):
            tbx, tbd1, tbd2, tbd3 = [], [], [], []
            for n in range(elementsCountAround + 1):
                tbx.append(nx[n])
                tbd1.append(nd1[n])
                tbd2.append([arcLengthAlong * vector.normalise(self._alongAxis)[c] for c in range(3)])
                tbd3.append(vector.normalise(vector.crossproduct3(tbd1[n], tbd2[n])))
            tnx.append(tbx)
            tnd1.append(tbd1)
            tnd2.append(tbd2)
            tnd3.append(tbd3)
            if self._cylinderType == CylinderType.CYLIDNER_TRUNCATED_CONE:
                nx = nd1 = []
                majorRadius1 = majorRadius1 * (1 - rate)
                majorAxis1 = vector.setMagnitude(self._majorAxis, majorRadius1)
                baseC = [self._baseCentre[c] + (n3 + 1) * arcLengthAlong * vector.normalise(self._alongAxis)[c] for c in
                         range(3)]
                nx, nd1 = self.createCylinderBaseMesh2D(
                    baseC, majorAxis1, minorAxis, elementsCountAround, majorRadius1)

        btx = self._shield.px
        btd1 = self._shield.pd1
        btd2 = self._shield.pd2
        btd3 = self._shield.pd3

        for n3 in range(self._elementsCountUp + 1):
            for n in range(self._shield.elementsCountAroundFull + 1):
                n1, n2 = self._shield.convertRimIndex(n)
                tx, td1, td2, td3 = tnx[n3], tnd1[n3], tnd2[n3], tnd3[n3]
                btx[n3][n2][n1] = tx[n]
                if n2 > self._shield.elementsCountRim:  # regular rows
                    btd1[n3][n2][n1] = td1[n]
                    btd3[n3][n2][n1] = td3[n]
                if n2 >= 2:
                    btd3[n3][n2][n1] = vector.setMagnitude(minorAxis, vector.dotproduct(td3[n], minorAxis))
                else:  # around rim
                    btd1[n3][n2][n1] = td1[n]
                    btd3[n3][n2][n1] = td3[n]
                btd2[n3][n2][n1] = td2[n]

        return btx, btd1, btd2, btd3

    def createMirrorCurve(self, n3):
        '''
        generate coordinates and derivatives for the mirror curve
        :param n3: Index of along cylinder axis coordinates to use
        :return: Coordinates and derivatives for the mirror curve
        '''
        btx = self._shield.px
        btd1 = self._shield.pd1
        btd2 = self._shield.pd2
        btd3 = self._shield.pd3

        rcx = []
        tmdx = btx[n3][0][(self._elementsCountAcross) // 2]
        tmdd3 = btd3[n3][0][(self._elementsCountAcross) // 2]
        tmux = [
            0.5 * (btx[n3][self._elementsCountUp][0][c] + btx[n3][self._elementsCountUp][self._elementsCountAcross][c])
            for c in range(3)]
        rcx.append(tmdx)
        rcx.append(tmux)
        rcd3 = [vector.setMagnitude(tmdd3, -1), vector.setMagnitude(tmdd3, -1)]
        rscx, rscd1 = sampleCubicHermiteCurves(rcx, rcd3, self._shield.elementsCountUp, lengthFractionStart=1,
                                               arcLengthDerivatives=True)[0:2]

        # get d2, d3
        rscd2 = []
        rscd3 = []
        for n in range(len(rscx)):
            d3 = vector.normalise(
                [btx[n3][self._elementsCountUp][self._elementsCountAcross][c] - btx[n3][self._elementsCountUp][0][c] for
                 c in range(3)])
            d2 = vector.normalise(vector.crossproduct3(d3, rscd1[n]))
            rscd2.append(d2)
            rscd3.append(d3)

        return rscx, rscd1, rscd2, rscd3

    def createRegularRowCurves(self, n3, rscx, rscd1, rscd3):
        '''
        generate curves along regular rows using the mirror curve obtained from createMirrorCurve.
        :param n3: Index of along cylinder axis coordinates to use
        '''
        btx = self._shield.px
        btd1 = self._shield.pd1
        btd2 = self._shield.pd2
        btd3 = self._shield.pd3

        elementsCountRim = 0
        for n2 in range(elementsCountRim + 2, self._elementsCountUp + 1):
            btx[n3][n2], btd3[n3][n2], pe, pxi, psf = sampleCubicHermiteCurves(
                [btx[n3][n2][0], rscx[n2], btx[n3][n2][-1]],
                [vector.setMagnitude(btd3[n3][n2][0], -1.0), rscd3[n2], btd3[n3][n2][-1]], self._elementsCountAcross,
                lengthFractionStart=1, lengthFractionEnd=1, arcLengthDerivatives=True)
            btd1[n3][n2] = \
            interpolateSampleCubicHermite([[-btd1[n3][n2][0][c] for c in range(3)], rscd1[n2], btd1[n3][n2][-1]],
                                          [[0.0, 0.0, 0.0]] * 3, pe, pxi, psf)[0]
            btd3[n3][n2][0] = [-btd3[n3][n2][0][c] for c in range(3)]
            btd1[n3][n2][0] = [-btd1[n3][n2][0][c] for c in range(3)]

    def createRegularColumnCurves(self, n3):
        '''
        up regular columns of shield: get d1, initial d3 below regular rows
        '''
        btx = self._shield.px
        btd1 = self._shield.pd1
        btd2 = self._shield.pd2
        btd3 = self._shield.pd3

        for n1 in range(2, self._elementsCountAcross - 1):
            tx, td1, pe, pxi, psf = sampleCubicHermiteCurves(
                [btx[n3][0][n1], btx[n3][2][n1]], [[-btd3[n3][0][n1][c] for c in range(3)], btd1[n3][2][n1]], 2,
                lengthFractionStart=1, arcLengthDerivatives=True)  # GRC fudge factor rvSulcusEdgeFactor
            for n2 in range(3, self._elementsCountUp + 1):
                tx.append(btx[n3][n2][n1])
                td1.append(btd1[n3][n2][n1])
            td1 = smoothCubicHermiteDerivativesLine(tx, td1, fixStartDirection=True)
            td3 = \
            interpolateSampleCubicHermite([btd1[n3][0][n1], btd3[n3][2][n1]], [[0.0, 0.0, 0.0]] * 2, pe, pxi, psf)[0]
            for n2 in range(self._elementsCountUp + 1):
                if n2 < 2:
                    btx[n3][n2][n1] = tx[n2]
                    if n2 == 0:
                        btd3[n3][n2][n1] = [-td1[0][c] for c in range(3)]
                    else:
                        btd3[n3][n2][n1] = td3[n2]
                if n2 == 0:
                    btd1[n3][n2][n1] = td3[n2]
                else:
                    btd1[n3][n2][n1] = td1[n2]

    def smoothTriplePointsCurves(self,n3,n2b,n1b,m1a):
        '''
        Smooth row and column curves passing triple points (i.e., row 1 and columns 1 and -2).
        '''
        btx = self._shield.px
        btd1 = self._shield.pd1
        btd2 = self._shield.pd2
        btd3 = self._shield.pd3

        # smooth shield row 1
        btd3[n3][n2b][n1b:m1a] = smoothCubicHermiteDerivativesLine(btx[n3][n2b][n1b:m1a], btd3[n3][n2b][n1b:m1a])

        # smooth Shield columns 1, -2
        for n1 in [1, -2]:
            tx = []
            td1 = []
            for n2 in range(1, self._elementsCountUp + 1):
                tx.append(btx[n3][n2][n1])
                td1.append(btd1[n3][n2][n1])
            td1 = smoothCubicHermiteDerivativesLine(tx, td1)
            for n in range(self._elementsCountUp):
                btd1[n3][n + 1][n1] = td1[n]

    def calculateD2Derivatives(self, n3, n3Count, arcLengthAlong):
        '''
        calculate d2 derivatives.
        :param n3: Index of along cylinder axis coordinates to use
        :param n3Count: number of bases to create coordinates for.
        :param arcLengthAlong: arc length along cylinder axis
        '''
        btx = self._shield.px
        btd1 = self._shield.pd1
        btd2 = self._shield.pd2
        btd3 = self._shield.pd3
        # get base d2 and next base x, d2
        for n2 in range(self._elementsCountUp + 1):
            for n1 in range(self._elementsCountAcross + 1):
                if btd1[n3][n2][n1]:
                    if n3 < n3Count:
                        btd2[n3][n2][n1] = [btx[n3 + 1][n2][n1][c] - btx[n3][n2][n1][c] for c in range(3)]
                    else:
                        btd2[n3][n2][n1] = vector.setMagnitude(vector.crossproduct3(btd3[n3][n2][n1], btd1[n3][n2][n1]),
                                                               arcLengthAlong)

    def smoothd2Derivatives(self):
        '''
        smooth d2 derivatives using initial values calculated by calculateD2Derivatives
        '''
        btx = self._shield.px
        btd1 = self._shield.pd1
        btd2 = self._shield.pd2
        btd3 = self._shield.pd3
        for n2 in range(self._elementsCountUp + 1):
            for n1 in range(self._elementsCountAcross + 1):
                td2 = []
                tx = []
                if btx[0][n2][n1]:
                    for n3 in range(self._elementsCountAlong + 1):
                        tx.append(btx[n3][n2][n1])
                        td2.append(btd2[n3][n2][n1])
                    td2 = smoothCubicHermiteDerivativesLine(tx, td2, fixStartDirection=True)
                    for n3 in range(self._elementsCountAlong + 1):
                        btd2[n3][n2][n1] = td2[n3]

    def createCylinderMesh3d(self, fieldModule, coordinates):
        """
        Create an extruded shape (ellipse/circle) mesh. Currently limited to ellipse or circle base with the alongAxis
        perpendicular to the base.
        :param fieldmodule: Zinc fieldmodule to create elements in.
        :param coordinates: Coordinate field to define.
        :return: Final values of nextNodeIdentifier, nextElementIdentifier.
        """
        assert (self._elementsCountAlong > 0), 'createCylinderMesh3d:  Invalid number of along elements'
        assert (self._elementsCountAcross > 4), 'createCylinderMesh3d: Invalid number of around elements'
        assert (self._elementsCountUp > 2), 'createCylinderMesh3d: Invalid number of up elements'
        assert (self._cylinderMode in [self._cylinderMode.CYLINDER_MODE_FULL, self._cylinderMode.CYLINDER_MODE_LOWER_HALF]), 'createCylinerMesh3d: Invalid cylinder mode.'
        plane=[-d for d in self._majorAxis]+[-vector.dotproduct(self._majorAxis,self._baseCentre)]

        nodes = fieldModule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh = fieldModule.findMeshByDimension(3)
        cache = fieldModule.createFieldcache()

        # create the base ellipse
        minorAxis = vector.setMagnitude(vector.crossproduct3(self._alongAxis, self._majorAxis), self._minorRadius)
        majorRadius = vector.magnitude(self._majorAxis)
        elementsCountAround = 2 * (self._elementsCountUp - 2) + self._elementsCountAcross

        # the bottom curve node coordinates and derivatives
        arcLengthAlong = vector.magnitude(self._alongAxis)/self._elementsCountAlong
        elementsCountRim = 0

        shieldMode = ShieldMode.SHIELD_MODE_FULL if self._cylinderMode is self._cylinderMode.CYLINDER_MODE_FULL else ShieldMode.SHIELD_MODE_LOWER_HALF
        self._shield = ShieldMesh(self._elementsCountAcross, self._elementsCountUp, elementsCountRim, None, self._elementsCountAlong, shieldMode)

        # generate bases mesh along cylinder axis
        btx, btd1, btd2, btd3 = self.generateBasesMesh(majorRadius, elementsCountAround, arcLengthAlong, minorAxis, self._radiusReductionRate)

        n3Count = 0 if self._cylinderType == CylinderType.CYLIDNER_REGULAR else self._elementsCountAlong
        for n3 in range(n3Count+1):
            rscx, rscd1, rscd2, rscd3 = self.createMirrorCurve(n3)
            self.createRegularRowCurves(n3, rscx, rscd1, rscd3)
            self.createRegularColumnCurves(n3)
            self._shield.getTriplePoints(n3)
            n1b = 1
            m1a = self._shield.elementsCountAcross
            m1b = m1a - 1
            m1c = m1a - 2
            n2b = 1
            self.smoothTriplePointsCurves(n3, n2b, n1b, m1a)
            self._shield.smoothDerivativesToTriplePoints(n3, fixAllDirections=True)

        for n3 in range(n3Count + 1):
            self.calculateD2Derivatives(n3, n3Count, arcLengthAlong)

        if self._cylinderType == CylinderType.CYLIDNER_TRUNCATED_CONE:
            self.smoothd2Derivatives()

        # The other bases.
        temx = []
        if self._cylinderType == CylinderType.CYLIDNER_REGULAR:
            for n2 in range(self._elementsCountUp + 1):
                for n3 in range(self._elementsCountAlong + 1):
                    for n1 in range(self._elementsCountAcross + 1):
                        if self._shield.px[0][n2][n1]:
                            temx = [self._shield.px[0][n2][n1][c] + n3*arcLengthAlong*vector.normalise(self._alongAxis)[c] for c in range(3)]
                            self._shield.px[n3][n2][n1]=temx
                            self._shield.pd1[n3][n2][n1]=self._shield.pd1[0][n2][n1]
                            self._shield.pd2[n3][n2][n1]=self._shield.pd2[0][n2][n1]
                            self._shield.pd3[n3][n2][n1]=self._shield.pd3[0][n2][n1]

        #################
        # Create nodes
        #################

        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        self._startNodeIdentifier = nodeIdentifier
        nodeIdentifier = self._shield.generateNodes(fieldModule, coordinates, nodeIdentifier, plane)
        self._endNodeIdentifier = nodeIdentifier

        #################
        # Create elements
        #################

        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)
        self._startElementIdentifier = elementIdentifier
        elementIdentifier = self._shield.generateElements(fieldModule, coordinates, elementIdentifier, [])
        self._endElementIdentifier = elementIdentifier




