"""
Utility functions for generating a solid cylinder (extruded ellipse/circle).
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

class CylinderMode(Enum):
    CYLINDER_MODE_FULL = 1        # full cylinder is createted
    CYLINDER_MODE_LOWER_HALF = 2  # describe

def createCylinderMesh3d(fieldModule, coordinates, baseCentre, alongAxis, majorAxis, minorRadius,
                         elementsCountAcross, elementsCountUp, elementsCountAlong, startNodeIdentifier, startElementIdentifier,
                         cylinderMode = CylinderMode.CYLINDER_MODE_FULL, useCrossDerivatives = False):
    """
    Create an extruded shape (ellipse/circle) mesh. Currently limited to ellipse or circle base with the alongAxis
    perpendicular to the base.
    :param fieldmodule: Zinc fieldmodule to create elements in.
    :param coordinates: Coordinate field to define.
    :param baseCentre: The centre of the base of the cylinder (ellipse/circle).
    :param alongAxis: The axis that the base is extruded along.
    :param majorAxis: The major axis of the base. Should be perpendicular to alongAxis
    :param minorRadius: The minor radius of the base.
    :param minorRadius: The minor radius of the base.
    :param elementsCountAcross: Number of elements across the base minor radius.
    :param elementsCountUp: Number of elements up central axis of shield. Must be at least 2 + elementsCountRim.
    :param cylinderMode: A value from enum CylinderMode specifying.
    :return: Final values of nextNodeIdentifier, nextElementIdentifier.
    """
    assert (elementsCountAlong > 0), 'createCylinderMesh3d:  Invalid number of along elements'
    assert (elementsCountAcross > 4), 'createCylinderMesh3d: Invalid number of around elements'
    assert (elementsCountUp > 2), 'createCylinderMesh3d: Invalid number of up elements'
    assert (cylinderMode in [cylinderMode.CYLINDER_MODE_FULL, cylinderMode.CYLINDER_MODE_LOWER_HALF]), 'createCylinerMesh3d: Invalid cylinder mode.'
    plane=[-d for d in majorAxis]+[-vector.dotproduct(majorAxis,baseCentre)]

    nodes = fieldModule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    mesh = fieldModule.findMeshByDimension(3)
    cache = fieldModule.createFieldcache()

    # create the base ellipse
    minorAxis = vector.setMagnitude(vector.crossproduct3(alongAxis, majorAxis), minorRadius)
    majorRadius = vector.magnitude(majorAxis)
    elementsCountAround = 2 * (elementsCountUp - 2) + elementsCountAcross
    nx, nd1 = createCylinderBaseMesh2D(
        baseCentre, majorAxis, minorAxis, elementsCountAround, majorRadius)

    # the bottom curve node coordinates and derivatives
    arcLengthAlong = vector.magnitude(alongAxis)/elementsCountAlong
    elementsCountRim = 0

    shieldMode = ShieldMode.SHIELD_MODE_FULL if cylinderMode is cylinderMode.CYLINDER_MODE_FULL else ShieldMode.SHIELD_MODE_LOWER_HALF
    tShield = ShieldMesh(elementsCountAcross, elementsCountUp, elementsCountRim, None, elementsCountAlong, shieldMode)

    tnx, tnd1, tnd2, tnd3 = [], [], [], []
    for n3 in range(elementsCountUp + 1):
        tbx, tbd1, tbd2, tbd3 = [], [], [], []
        for n in range(elementsCountAround + 1):
            tbx.append(nx[n])
            tbd1.append(nd1[n])
            tbd2.append([arcLengthAlong * vector.normalise(alongAxis) [c] for c in range(3)])
            tbd3.append(vector.normalise(vector.crossproduct3(tbd1[n], tbd2[n])))
        tnx.append(tbx)
        tnd1.append(tbd1)
        tnd2.append(tbd2)
        tnd3.append(tbd3)

    btx  = tShield.px
    btd1 = tShield.pd1
    btd2 = tShield.pd2
    btd3 = tShield.pd3

    for n3 in range(elementsCountUp+1):
        for n in range(tShield.elementsCountAroundFull+1):
            n1, n2 = tShield.convertRimIndex(n)
            tx, td1, td2, td3 = tnx[n3], tnd1[n3], tnd2[n3], tnd3[n3]
            btx[n3][n2][n1] = tx[n]
            if n2 > tShield.elementsCountRim:  # regular rows
                # if n1 < tShield.elementsCountAcross:
                btd1[n3][n2][n1] = td1[n]
                btd3[n3][n2][n1] = td3[n]
            if n2 >= 2:
                btd3[n3][n2][n1] = vector.setMagnitude(minorAxis, vector.dotproduct(td3[n], minorAxis))
            else:  # around rim
                btd1[n3][n2][n1] = td1[n]
                btd3[n3][n2][n1] = td3[n]
            btd2[n3][n2][n1] = td2[n]

    # create the curve on the major axis
    rcx = []
    tmdx = btx[0][0][(elementsCountAcross)//2]
    tmdd3 = btd3[0][0][(elementsCountAcross)//2]
    tmux = [0.5*(btx[0][elementsCountUp][0][c]+btx[0][elementsCountUp][elementsCountAcross][c]) for c in range(3)]
    rcx.append(tmdx)
    rcx.append(tmux)
    rcd3 = [vector.setMagnitude(tmdd3, -1), vector.setMagnitude(tmdd3, -1)]
    rscx, rscd1 = sampleCubicHermiteCurves(rcx, rcd3, tShield.elementsCountUp, lengthFractionStart=0.667, arcLengthDerivatives = True)[0:2]

    # get d2, d3
    rscd2 = []
    rscd3 = []
    for n in range(len(rscx)):
        d3 = vector.normalise([tbx[tShield.elementsCountAroundFull][c] - tbx[0][c] for c in range(3)])
        d2 = vector.normalise(vector.crossproduct3(d3, rscd1[n]))
        rscd2.append(d2)
        rscd3.append(d3)

    # across regular rows of Shield: get d3, initial d1.
    for n2 in range(elementsCountRim + 2, elementsCountUp + 1):
        btx[0][n2], btd3[0][n2], pe, pxi, psf = sampleCubicHermiteCurves(
            [ btx[0][n2][0], rscx[n2], btx[0][n2][-1] ], [ vector.setMagnitude(btd3[0][n2][0], -1.0), rscd3[n2], btd3[0][n2][-1] ], elementsCountAcross,
            lengthFractionStart=0.667, lengthFractionEnd=0.667, arcLengthDerivatives = True)
        btd1[0][n2] = interpolateSampleCubicHermite([[-btd1[0][n2][0][c] for c in range(3)], rscd1[n2], btd1[0][n2][-1]], [[0.0, 0.0, 0.0]] * 3, pe, pxi, psf)[0]
        btd3[0][n2][0] = [-btd3[0][n2][0][c] for c in range(3)]
        btd1[0][n2][0] = [-btd1[0][n2][0][c] for c in range(3)]

    # up regular columns of shield: get d1, initial d3 below regular rows
    for n1 in range(2, elementsCountAcross - 1):
        tx, td1, pe, pxi, psf = sampleCubicHermiteCurves(
            [ btx[0][0][n1], btx[0][2][n1] ], [ [-btd3[0][0][n1][c] for c in range(3)], btd1[0][2][n1] ], 2, lengthFractionStart=0.667, arcLengthDerivatives = True)  # GRC fudge factor rvSulcusEdgeFactor
        for n2 in range(3, elementsCountUp + 1):
            tx .append(btx [0][n2][n1])
            td1.append(btd1[0][n2][n1])
        td1 = smoothCubicHermiteDerivativesLine(tx, td1, fixStartDirection = True)
        td3 = interpolateSampleCubicHermite([ btd1[0][0][n1], btd3[0][2][n1] ], [ [ 0.0, 0.0, 0.0 ] ]*2, pe, pxi, psf)[0]
        for n2 in range(elementsCountUp + 1):
            if n2 < 2:
                btx [0][n2][n1] = tx [n2]
                if n2 == 0:
                    btd3[0][n2][n1] = [-td1[0][c] for c in range(3) ]
                else:
                    btd3[0][n2][n1] = td3[n2]
            if n2 == 0:
                btd1[0][n2][n1] = td3[n2]
            else:
                btd1[0][n2][n1] = td1[n2]

    tShield.getTriplePoints(n3=0)
    n1b = 1
    m1a = tShield.elementsCountAcross
    m1b = m1a - 1
    m1c = m1a - 2
    n2b = 1

    # smooth shield row 1
    btd3[0][n2b][n1b:m1a] = smoothCubicHermiteDerivativesLine(btx[0][n2b][n1b:m1a], btd3[0][n2b][n1b:m1a])

    # smooth Shield columns 1, -2
    for n1 in [ 1, -2 ]:
        tx = []
        td1 = []
        for n2 in range(1, elementsCountUp + 1):
            tx .append(btx [0][n2][n1])
            td1.append(btd1[0][n2][n1])
        td1 = smoothCubicHermiteDerivativesLine(tx, td1)
        for n in range(elementsCountUp):
            btd1[0][n + 1][n1] = td1[n]

    tShield.smoothDerivativesToTriplePoints(n3=0, fixAllDirections=True)

    # get base d2 and next base x, d2
    for n2 in range(elementsCountUp + 1):
        for n1 in range(elementsCountAcross + 1):
            if btd1[0][n2][n1]:
                btd2[0][n2][n1] = btd2[1][n2][n1] = vector.setMagnitude(vector.crossproduct3(btd3[0][n2][n1], btd1[0][n2][n1]), arcLengthAlong)
                btx [1][n2][n1] = [ (btx [0][n2][n1][c] + btd2[0][n2][n1][c]) for c in range(3) ]

    # get next base d1, d3
    # row 1
    btd3[1][n2b][n1b:m1a] = smoothCubicHermiteDerivativesLine(btx[1][n2b][n1b:m1a], btd3[0][n2b][n1b:m1a], fixAllDirections = True)
    # regular rows 2+
    for n2 in range(2, elementsCountUp + 1):
        btd3[1][n2] = smoothCubicHermiteDerivativesLine(btx[1][n2], [[-d for d in btd3[0][n2][0]]]+btd3[0][n2][1:], fixAllDirections = True)
        btd3[1][n2][0] = [-d for d in btd3[1][n2][0]]
    # columns
    for n1 in range(n1b, m1a):
        startn2 = 1 if (n1 in [n1b, m1b]) else 0
        tx  = []
        td1 = []
        for n2 in range(startn2, elementsCountUp + 1):
            tx .append(btx [1][n2][n1])
            if n2 == 0:
                td1.append([-d for d in btd3[0][n2][n1]])
            else:
                td1.append(btd1[0][n2][n1])
        td1 = smoothCubicHermiteDerivativesLine(tx, td1, fixAllDirections = True)
        for n2 in range(startn2, elementsCountUp + 1):
            if n2 == 0:
                btd3[1][n2][n1] = [-d for d in td1[0]]
            else:
                btd1[1][n2][n1] = td1[n2 - startn2]

    # fix next base derivatives leading to triple points
    # first copy d1 from base to next base
    for n1 in [ n1b, m1b ]:
        for n2 in range(0, n2b):
            btd3[1][n2][n1] = btd3[0][n2][n1]
    tShield.smoothDerivativesToTriplePoints(n3=1, fixAllDirections=True)

    # btd3[1][0][1] = smoothCubicHermiteDerivativesLine([ btx[1][0][1], btx[1][1][1] ], [ [-d for d in btd3[1][0][1]], [ (btd1[1][1][1][c] + btd2[1][1][1][c]) for c in range(3) ] ],
    #                                                   fixEndDerivative = True, fixStartDirection = True)[0]
    # btd2[1][0][m1b] = smoothCubicHermiteDerivativesLine([ btx[1][0][m1b], btx[1][1][m1b] ], [ btd2[1][0][m1b], [ (-btd1[1][1][m1b][c] + btd2[1][1][m1b][c]) for c in range(3) ] ],
    #                                                   fixEndDerivative = True, fixStartDirection = True)[0]

    # The other bases.
    temx = []
    if elementsCountAlong>1:
        for n2 in range(elementsCountUp + 1):
            for n3 in range(2, elementsCountAlong + 1):
                for n1 in range(elementsCountAcross + 1):
                    if tShield.px[0][n2][n1]:
                        temx = [tShield.px[0][n2][n1][c] + n3*arcLengthAlong*vector.normalise(alongAxis)[c] for c in range(3)]
                        tShield.px[n3][n2][n1]=temx
                        tShield.pd1[n3][n2][n1]=tShield.pd1[0][n2][n1]
                        tShield.pd2[n3][n2][n1]=tShield.pd2[0][n2][n1]
                        tShield.pd3[n3][n2][n1]=tShield.pd3[0][n2][n1]

    #################
    # Create nodes
    #################

    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

    nodeIdentifier = startNodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
    nodeIdentifier = tShield.generateNodes(fieldModule, coordinates, nodeIdentifier, plane)

    #################
    # Create elements
    #################

    elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)
    elementIdentifier = tShield.generateElements(fieldModule, coordinates, elementIdentifier, [])




