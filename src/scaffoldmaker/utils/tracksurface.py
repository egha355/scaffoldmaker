'''
Utility class for representing surfaces on which features can be located.
'''

from __future__ import division
import copy
from enum import Enum
import math
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import vector


class TrackSurface:
    '''
    A surface description on which positions can be stored and tracked for
    specified directions and distances for location surface features.
    Currently represented by a lattice of elementsCount1 x elementsCount2
    square elements with bicubic Hermite interpolation but zero cross derivatives.
    '''

    def __init__(self, elementsCount1, elementsCount2, nx, nd1, nd2, loop1 = False):
        '''
        Creates a TrackSurface with a lattice of elementsCount1*elementsCount2
        elements, with nodes and elements varying across direction 1 fastest.
        :param elementsCount1, elementsCount2: Number of elements in directions 1 and 2.
        :param nx, nd1, nd2: List of (elementsCount2 + 1)*(elementsCount1 + 1) node
        coordinates, derivatives 1 and derivatives 2. Cross derivatives are zero.
        All coordinates are 3 component.
        :param loop1: Set to True if loops back to start in direction 1. Means:
        - one fewer node supplied around
        - all nodes added twice plus first node a third time
        - twice as many elements
        - valid proportions are from 0.0 to 2.0.
        '''
        self.elementsCount1 = elementsCount1
        self.elementsCount2 = elementsCount2
        self.nx  = []
        self.nd1 = []
        self.nd2 = []
        if loop1:
            n = 0
            for n2 in range(elementsCount2 + 1):
                nLimit = n + elementsCount1
                rx  = nx [n:nLimit]
                rd1 = nd1[n:nLimit]
                rd2 = nd2[n:nLimit]
                self.nx  += rx  + rx  + [ rx [0] ]
                self.nd1 += rd1 + rd1 + [ rd1[0] ]
                self.nd2 += rd2 + rd2 + [ rd2[0] ]
                n += elementsCount1
            self.elementsCount1 *= 2
        else:
            self.nx  += nx
            self.nd1 += nd1
            self.nd2 += nd2
        self.loop1 = loop1

    def createMirrorX(self):
        '''
        Mirror track surface about x axis by negating all x coordinates
        and rewind elements by flipping order and derivatives in direction 1.
        :return: TrackSurface
        '''
        nx  = []
        nd1 = []
        nd2 = []
        nodesCount1 = self.elementsCount1 + (0 if self.loop1 else 1)
        nodesCount2 = self.elementsCount2 + 1
        for n2 in range(nodesCount2):
            for n1 in range(nodesCount1):
                oi = n2*nodesCount1 + self.elementsCount1 - n1
                ox  = copy.deepcopy(self.nx [oi])
                od1 = copy.deepcopy(self.nd1[oi])
                od2 = copy.deepcopy(self.nd2[oi])
                nx .append([ -ox [0],  ox [1],  ox [2] ])
                nd1.append([  od1[0], -od1[1], -od1[2] ])
                nd2.append([ -od2[0],  od2[1],  od2[2] ])
        return TrackSurface(self.elementsCount1, self.elementsCount2, nx, nd1, nd2, loop1 = self.loop1)

    def createPositionProportion(self, proportion1, proportion2):
        '''
        Return position on surface for proportions across directions 1 and 2.
        :param proportion1, proportion2: Proportions across directions 1 and 2,
        each varying from 0.0 to 1.0 (or 2.0 if loop), with elements equal sized.
        :return: TrackSurfacePosition
        '''
        maxProportion1 = 2.0 if self.loop1 else 1.0
        assert (proportion1 >= 0.0) and (proportion1 <= maxProportion1), 'createPositionProportion:  Proportion 1 (' + str(proportion1) + ') out of range'
        assert (proportion2 >= 0.0) and (proportion2 <= 1.0), 'createPositionProportion:  Proportion 2 (' + str(proportion2) + ') out of range'

        pe1 = (proportion1/maxProportion1)*self.elementsCount1
        if pe1 < self.elementsCount1:
            e1 = int(pe1)
            xi1 = pe1 - e1
        else:
            e1 = self.elementsCount1 - 1
            xi1 = 1.0
        pe2 = proportion2*self.elementsCount2
        if pe2 < self.elementsCount2:
            e2 = int(pe2)
            xi2 = pe2 - e2
        else:
            e2 = self.elementsCount2 - 1
            xi2 = 1.0
        return TrackSurfacePosition(e1, e2, xi1, xi2)

    def getProportion(self, position):
        '''
        From a position on this track surface, return proportions.
        :return: proportion1, proportion2
        '''
        maxProportion1 = 2.0 if self.loop1 else 1.0
        return [ maxProportion1*(position.e1 + position.xi1)/self.elementsCount1, (position.e2 + position.xi2)/self.elementsCount2 ]

    def evaluateCoordinates(self, position, derivatives = False):
        '''
        Evaluate coordinates on surface at position, and optionally
        derivative w.r.t. xi1 and xi2.
        :param position: A valid TrackSurfacePosition.
        :return: If derivatives is False: coordinates [ x, y, z].
        If derivatives is True: coordinates, derivative1, derivative2.
        '''
        n1 = position.e2*(self.elementsCount1 + 1) + position.e1
        n2 = n1 + 1
        n3 = n1 + self.elementsCount1 + 1
        n4 = n3 + 1
        ni = [ n1, n2, n3, n4 ]
        f1x1, f1d1, f1x2, f1d2 = interp.getCubicHermiteBasis(position.xi1)
        f2x1, f2d1, f2x2, f2d2 = interp.getCubicHermiteBasis(position.xi2)
        fx  = [ f1x1*f2x1, f1x2*f2x1, f1x1*f2x2, f1x2*f2x2 ]
        fd1 = [ f1d1*f2x1, f1d2*f2x1, f1d1*f2x2, f1d2*f2x2 ]
        fd2 = [ f1x1*f2d1, f1x2*f2d1, f1x1*f2d2, f1x2*f2d2 ]
        nx = self.nx
        nd1 = self.nd1
        nd2 = self.nd2
        coordinates = []
        for c in range(3):
            x = 0.0
            for ln in range(4):
                gn = ni[ln]
                x += fx[ln]*nx[gn][c] + fd1[ln]*nd1[gn][c] + fd2[ln]*nd2[gn][c]
            coordinates.append(x)
        if not derivatives:
            return coordinates
        df1x1, df1d1, df1x2, df1d2 = interp.getCubicHermiteBasisDerivatives(position.xi1)
        d1fx  = [ df1x1*f2x1, df1x2*f2x1, df1x1*f2x2, df1x2*f2x2 ]
        d1fd1 = [ df1d1*f2x1, df1d2*f2x1, df1d1*f2x2, df1d2*f2x2 ]
        d1fd2 = [ df1x1*f2d1, df1x2*f2d1, df1x1*f2d2, df1x2*f2d2 ]
        df2x1, df2d1, df2x2, df2d2 = interp.getCubicHermiteBasisDerivatives(position.xi2)
        d2fx  = [ f1x1*df2x1, f1x2*df2x1, f1x1*df2x2, f1x2*df2x2 ]
        d2fd1 = [ f1d1*df2x1, f1d2*df2x1, f1d1*df2x2, f1d2*df2x2 ]
        d2fd2 = [ f1x1*df2d1, f1x2*df2d1, f1x1*df2d2, f1x2*df2d2 ]
        derivative1 = []
        derivative2 = []
        for c in range(3):
            d1 = 0.0
            d2 = 0.0
            for ln in range(4):
                gn = ni[ln]
                d1 += d1fx[ln]*nx[gn][c] + d1fd1[ln]*nd1[gn][c] + d1fd2[ln]*nd2[gn][c]
                d2 += d2fx[ln]*nx[gn][c] + d2fd1[ln]*nd1[gn][c] + d2fd2[ln]*nd2[gn][c]
            derivative1.append(d1)
            derivative2.append(d2)
        return coordinates, derivative1, derivative2

    class HermiteCurveMode(Enum):
        SMOOTH = 1    # smooth variation of element size between end derivatives
        TRANSITION_END = 2  # transition from start derivative then even size
        TRANSITION_START = 3  # even size then transition to end derivative
        TRANSITION_START_AND_END = 3  # transition to/from start and end derivatives and even size in between
        UNIFORM_SIZE = 4  # regardless of initial derivatives, sample in even sizes

    def createHermiteCurvePoints(self, aProportion1, aProportion2, bProportion1, bProportion2, elementsCount,
            derivativeStart = None, derivativeEnd = None, curveMode = HermiteCurveMode.SMOOTH):
        '''
        Create hermite curve points between two points a and b on the surface, each defined
        by their proportions over the surface in directions 1 and 2.
        Also returns cross direction 2 in plane of surface with similar magnitude to curve derivative 1,
        and unit surface normals.
        :param derivativeStart, derivativeEnd: Optional derivative vectors in 3-D world coordinates
        to match at the start and end of the curves. If omitted, fits in with other derivative or is
        in a straight line from a to b.
        :param elementsCount:  Number of elements out.
        :return: nx[], nd1[], nd2[], nd3[], nProportions[]
        '''
        #print('createHermiteCurvePoints', aProportion1, aProportion2, bProportion1, bProportion2, elementsCount, derivativeStart, derivativeEnd)
        if derivativeStart:
            position = self.createPositionProportion(aProportion1, aProportion2)
            _, sd1, sd2 = self.evaluateCoordinates(position, derivatives = True)
            delta_xi1, delta_xi2 = calculate_surface_delta_xi(sd1, sd2, derivativeStart)
            dp1Start = delta_xi1/self.elementsCount1
            if self.loop1:
                dp1Start *= 2.0
            dp2Start = delta_xi2/self.elementsCount2
            derivativeMagnitudeStart = math.sqrt(dp1Start*dp1Start + dp2Start*dp2Start)
            dp1Start *= elementsCount
            dp2Start *= elementsCount
            #print('start delta_xi1', delta_xi1, 'delta_xi2', delta_xi2)
            #print('dp1Start', dp1Start, 'dp2Start', dp2Start)
        if derivativeEnd:
            position = self.createPositionProportion(bProportion1, bProportion2)
            _, sd1, sd2 = self.evaluateCoordinates(position, derivatives = True)
            delta_xi1, delta_xi2 = calculate_surface_delta_xi(sd1, sd2, derivativeEnd)
            dp1End = delta_xi1/self.elementsCount1
            dp2End = delta_xi2/self.elementsCount2
            derivativeMagnitudeEnd = math.sqrt(dp1End*dp1End + dp2End*dp2End)
            dp1End *= elementsCount
            if self.loop1:
                dp1End *= 2.0
            dp2End *= elementsCount
            #print('end delta_xi1', delta_xi1, 'delta_xi2', delta_xi2)
            #print('dp1End', dp1End, 'dp2End', dp2End)
        if not derivativeStart:
            if derivativeEnd:
                dp1Start, dp2Start = interp.interpolateLagrangeHermiteDerivative([ aProportion1, aProportion2 ], [ bProportion1, bProportion2 ], [ dp1End, dp2End ], 0.0)
            else:
                dp1Start = bProportion1 - aProportion1
                dp2Start = bProportion2 - aProportion2
            derivativeMagnitudeStart = math.sqrt(dp1Start*dp1Start + dp2Start*dp2Start)/elementsCount
        if not derivativeEnd:
            if derivativeStart:
                dp1End, dp2End = interp.interpolateHermiteLagrangeDerivative([ aProportion1, aProportion2 ], [ dp1Start, dp2Start ], [ bProportion1, bProportion2 ], 1.0)
            else:
                dp1End = bProportion1 - aProportion1
                dp2End = bProportion2 - aProportion2
            derivativeMagnitudeEnd = math.sqrt(dp1End*dp1End + dp2End*dp2End)/elementsCount
        maxProportion1 = 2.0 if self.loop1 else 1.0
        #print('derivativeMagnitudeStart', derivativeMagnitudeStart, 'derivativeMagnitudeEnd', derivativeMagnitudeEnd)
        proportions, dproportions = interp.sampleCubicHermiteCurvesSmooth([ [ aProportion1, aProportion2 ], [ bProportion1, bProportion2 ] ], \
            [ [ dp1Start, dp2Start ], [ dp1End, dp2End ] ], elementsCount, derivativeMagnitudeStart, derivativeMagnitudeEnd)[0:2]
        if curveMode != self.HermiteCurveMode.SMOOTH:
            if derivativeStart and (curveMode in [self.HermiteCurveMode.TRANSITION_START, self.HermiteCurveMode.TRANSITION_START_AND_END]):
                addLengthStart = 0.5*derivativeMagnitudeStart
                lengthFractionStart = 0.5
            else:
                addLengthStart = 0.0
                lengthFractionStart = 1.0
            if derivativeEnd and (curveMode in [self.HermiteCurveMode.TRANSITION_END, self.HermiteCurveMode.TRANSITION_START_AND_END]):
                addLengthEnd = 0.5*derivativeMagnitudeEnd
                lengthFractionEnd = 0.5
            else:
                addLengthEnd = 0.0
                lengthFractionEnd = 1.0
            proportions, dproportions = interp.sampleCubicHermiteCurves(proportions, dproportions, elementsCount,
                                                                        addLengthStart, addLengthEnd, lengthFractionStart, lengthFractionEnd)[0:2]
        #print(' proportions', proportions)
        #print('dproportions', dproportions)
        nx  = []
        nd1 = []
        nd2 = []
        nd3 = []
        for n in range(0, elementsCount + 1):
            position = self.createPositionProportion(proportions[n][0], proportions[n][1])
            x, sd1, sd2 = self.evaluateCoordinates(position, derivatives = True)
            f1 = dproportions[n][0]*self.elementsCount1/maxProportion1
            f2 = dproportions[n][1]*self.elementsCount2
            d1 = [ (f1*sd1[c] + f2*sd2[c]) for c in range(3) ]
            d3 = vector.crossproduct3(sd1, sd2)
            # handle zero magnitude of d3
            mag = math.sqrt(sum(d3[c]*d3[c] for c in range(3)))
            if mag > 0.0:
                d3 = [ (d3[c]/mag) for c in range(3) ]
            d2 = vector.crossproduct3(d3, d1)
            nx .append(x)
            nd2.append(d2)
            nd1.append(d1)
            nd3.append(d3)
        #print('createHermiteCurvePoints end \n nx', nx,'\nnd1',nd1,'\nnd2',nd2,'\nnd3',nd3)
        return nx, nd1, nd2, nd3, proportions

    def resampleHermiteCurvePointsSmooth(self, nx, nd1, nd2, nd3, nProportions, derivativeMagnitudeStart=None, derivativeMagnitudeEnd=None):
        '''
        Call interp.sampleCubicHermiteCurvesSmooth on nx, nd1 and recalculate positions, nd2, nd3 for points.
        :return: nx[], nd1[], nd2[], nd3[], nProportions[]
        '''
        elementsCount = len(nx) - 1
        #print(nx, nd1, elementsCount, derivativeMagnitudeStart, derivativeMagnitudeEnd)
        nx, nd1 = interp.sampleCubicHermiteCurvesSmooth(nx, nd1, elementsCount, derivativeMagnitudeStart, derivativeMagnitudeEnd)[0:2]
        mag2 = vector.magnitude(nd2[0])
        if mag2 > 0.0:
            nd2[0] = vector.setMagnitude(nd2[0], vector.magnitude(nd1[0]))
        for n in range(1, elementsCount):
            p = self.findNearestPosition(nx[n], self.createPositionProportion(*nProportions[n]))
            nProportions[n] = self.getProportion(p)
            _, sd1, sd2 = self.evaluateCoordinates(p, derivatives=True)
            _, d2, d3 = calculate_surface_axes(sd1, sd2, vector.normalise(nd1[n]))
            nd2[n] = vector.setMagnitude(d2, vector.magnitude(nd1[n]))
            nd3[n] = d3
        mag2 = vector.magnitude(nd2[-1])
        if mag2 > 0.0:
            nd2[-1] = vector.setMagnitude(nd2[-1], vector.magnitude(nd1[-1]))
        return nx, nd1, nd2, nd3, nProportions

    def findNearestPosition(self, targetx, startPosition = None):
        '''
        Find the nearest point to targetx on the track surface, with optional start position.
        Only works if track surface is simply shaped; use a close startPosition if not.
        :return: Nearest TrackSurfacePosition
        '''
        if not startPosition:
            startPosition = self.createPositionProportion(0.5, 0.5)
        position = copy.deepcopy(startPosition)
        max_mag_dxi = 0.5  # target/maximum magnitude of xi increment
        xi_tol = 1.0E-6
        for iter in range(100):
            xi1 = position.xi1
            xi2 = position.xi2
            ax, ad1, ad2 = self.evaluateCoordinates(position, derivatives = True)
            deltax = [ (targetx[c] - ax[c]) for c in range(3) ]
            #print('iter', iter + 1, 'position', position, 'deltax', deltax, 'err', vector.magnitude(deltax))
            adelta_xi1, adelta_xi2 = calculate_surface_delta_xi(ad1, ad2, deltax)
            mag_dxi = math.sqrt(adelta_xi1*adelta_xi1 + adelta_xi2*adelta_xi2)
            dxi1 = adelta_xi1
            dxi2 = adelta_xi2
            if mag_dxi > max_mag_dxi:
                dxi1 *= max_mag_dxi/mag_dxi
                dxi2 *= max_mag_dxi/mag_dxi
            #print('    dxi', dxi1, dxi2)
            bxi1, bxi2, proportion, faceNumber = increment_xi_on_square(xi1, xi2, dxi1, dxi2)
            position.xi1 = bxi1
            position.xi2 = bxi2
            if mag_dxi < xi_tol:
                #print('converged mag_dxi', mag_dxi)
                break
            if faceNumber:
                onBoundary = self.updatePositionTofaceNumber(position, faceNumber)
                if onBoundary and (proportion < xi_tol):
                    # slide along boundary to nearest point; may cross other sides
                    if faceNumber in [1, 2]:
                        bdxi1 = 0.0
                        cmag_dxi = bdxi2 = dxi2*(1.0 - proportion)
                    else:
                        cmag_dxi = bdxi1 = dxi1*(1.0 - proportion)
                        bdxi2 = 0.0
                    cxi1, cxi2, cproportion, cFaceNumber = increment_xi_on_square(bxi1, bxi2, bdxi1, bdxi2)
                    position.xi1 = cxi1
                    position.xi2 = cxi2
                    if math.fabs(cmag_dxi) < xi_tol:
                        #print('converged on boundary cmag_dxi', cmag_dxi)
                        break
                    if cFaceNumber:
                        cOnBoundary = self.updatePositionTofaceNumber(position, cFaceNumber)
                        if cOnBoundary:
                            #print('TrackSurface.findNearestPosition:  End on corner boundary at', position)
                            break
        else:
            print('TrackSurface.findNearestPosition:  Reach max iterations', iter + 1, 'closeness in xi', mag_dxi)
        #print('final position', position)
        return position

    def trackVector(self, startPosition, direction, trackDistance):
        '''
        Track from startPosition the given distance in the vector direction.
        Approximate, uses improved Euler method (mean of original & final gradient).
        :param startPosition: TrackSurfacePosition
        :param direction: 3-D vector (x, y, z) to track along. Projected onto surface.
        :param trackDistance: Distance to track along. Can be negative.
        :return: Final TrackSurfacePosition
        '''
        #print('TrackSurface.trackVector  start position', startPosition, 'direction', direction, 'distance', trackDistance)
        useDirection = direction
        useTrackDistance = trackDistance
        if trackDistance < 0.0:
            useDirection = [ -d for d in direction ]
            useTrackDistance = -trackDistance
        position = copy.deepcopy(startPosition)
        distance = 0.0
        distanceLimit = 0.9999*useTrackDistance
        max_mag_dxi = 0.02  # target/maximum magnitude of xi increment

        while distance < useTrackDistance:
            xi1 = position.xi1
            xi2 = position.xi2
            ax, ad1, ad2 = self.evaluateCoordinates(position, derivatives = True)
            adelta_xi1, adelta_xi2 = calculate_surface_delta_xi(ad1, ad2, useDirection)
            #print('adelta_xi', adelta_xi1, adelta_xi2)
            scale = max_mag_dxi/math.sqrt(adelta_xi1*adelta_xi1 + adelta_xi2*adelta_xi2)
            adxi1 = dxi1 = scale*adelta_xi1
            adxi2 = dxi2 = scale*adelta_xi2
            #print('adxi', adxi1, adxi2)
            for i in range(1):  # can try more in future
                # can go slightly outside element to get predictor/correct position
                position.xi1 = xi1 + dxi1
                position.xi2 = xi2 + dxi2
                bx, bd1, bd2 = self.evaluateCoordinates(position, derivatives = True)
                bdelta_xi1, bdelta_xi2 = calculate_surface_delta_xi(bd1, bd2, useDirection)
                # use mean of start and end derivatives
                delta_xi1 = 0.5*(adelta_xi1 + bdelta_xi1)
                delta_xi2 = 0.5*(adelta_xi2 + bdelta_xi2)
                scale = max_mag_dxi/math.sqrt(delta_xi1*delta_xi1 + delta_xi2*delta_xi2)
                dxi1 = scale*delta_xi1
                dxi2 = scale*delta_xi2
            bxi1, bxi2, proportion, faceNumber = increment_xi_on_square(xi1, xi2, dxi1, dxi2)
            position.xi1 = bxi1
            position.xi2 = bxi2
            #print(distance, '-->', position)
            bx, bd1, bd2 = self.evaluateCoordinates(position, derivatives = True)
            bdelta_xi1, bdelta_xi2 = calculate_surface_delta_xi(bd1, bd2, useDirection)
            scale = max_mag_dxi/math.sqrt(bdelta_xi1*bdelta_xi1 + bdelta_xi2*bdelta_xi2)
            bdxi1 = scale*bdelta_xi1
            bdxi2 = scale*bdelta_xi2
            #print('bdxi', bdxi1, bdxi2)
            ad = [ proportion*(adxi1*ad1[c] + adxi2*ad2[c]) for c in range(3) ]
            bd = [ proportion*(bdxi1*bd1[c] + bdxi2*bd2[c]) for c in range(3) ]
            # GRC check:
            #arcLength = interp.computeCubicHermiteArcLength(ax, ad, bx, bd, rescaleDerivatives = True)
            arcLength = interp.getCubicHermiteArcLength(ax, ad, bx, bd)
            #print('scales', vector.magnitude([ (bx[c] - ax[c]) for c in range(3) ]), vector.magnitude(ad), vector.magnitude(bd), 'arc length', arcLength)
            if (distance + arcLength) >= distanceLimit:
                # limit to useTrackDistance, approximately, and finish
                r = proportion*(useTrackDistance - distance)/arcLength
                position.xi1 = xi1 + r*dxi1
                position.xi2 = xi2 + r*dxi2
                #print(distance, '-->', position, '(final)')
                break
            if (arcLength == 0.0) and (not faceNumber):
                print('TrackSurface.trackVector. No increment at', position, 'final distance', distance, 'of', useTrackDistance)
                break
            distance += arcLength
            if faceNumber:
                onBoundary = self.updatePositionTofaceNumber(position, faceNumber)
                if onBoundary:
                    print('TrackSurface.trackVector:  End on boundary at', position)
                    break
                #print('  cross face', faceNumber, 'new position', position)
        return position

    def updatePositionTofaceNumber(self, position, faceNumber):
        '''
        Update coordinates of TrackSurfacePosition position to cross
        the given face number, or clamp to range if reached boundary.
        :return: True if reached boundary of track surface, otherwise False.
        '''
        onBoundary = False
        if faceNumber == 1:  # xi1 == 0.0
            if position.e1 > 0:
                position.e1 -= 1
                position.xi1 = 1.0
            else:
                position.xi1 = 0.0
                onBoundary = True
        elif faceNumber == 2:  # xi1 == 1.0
            if position.e1 < (self.elementsCount1 - 1):
                position.e1 += 1
                position.xi1 = 0.0
            else:
                position.xi1 = 1.0
                onBoundary = True
        elif faceNumber == 3:  # xi2 == 0.0
            if position.e2 > 0:
                position.e2 -= 1
                position.xi2 = 1.0
            else:
                position.xi2 = 0.0
                onBoundary = True
        elif faceNumber == 4:  # xi2 == 1.0
            if position.e2 < (self.elementsCount2 - 1):
                position.e2 += 1
                position.xi2 = 0.0
            else:
                position.xi2 = 1.0
                onBoundary = True
        #if onBoundary:
        #    print('!!! Reached boundary of face', faceNumber, 'position', position)
        return onBoundary


class TrackSurfacePosition:
    '''
    Position on TrackSurface. Create with createPosition~ method of TrackSurface.
    '''

    def __init__(self, e1, e2, xi1, xi2):
        '''
        :param e1, e2: Element index in directions 1 and 2, starting at 0.
        :param xi1, xi2: 2-D element chart coordinates 0.0 <= xiN <= 1.0
        '''
        self.e1 = e1
        self.e2 = e2
        self.xi1 = xi1
        self.xi2 = xi2

    def __str__(self):
        return 'element (' + str(self.e1) + ',' + str(self.e2) + ') xi (' + str(self.xi1) + ',' + str(self.xi2) + ')'

def calculate_surface_delta_xi(d1, d2, direction):
    '''
    Calculate dxi1, dxi2 in 3-D vector direction.
    :param d1, d2: Derivatives of coordinate w.r.t. xi1, xi2.
    :param direction: 3-D vector.
    :return: delta_xi1, delta_xi2
    '''
    # overdetermined system (3-D vector, 2-D surface), solve least squares
    # A transpose A x = A transpose b
    a = [ [ 0.0, 0.0 ], [ 0.0, 0.0 ]]
    b = [ 0.0, 0.0 ]
    dx_dxi = [ d1, d2 ]
    for i in range(2):
        for j in range(2):
            for k in range(3):
                a[i][j] += dx_dxi[i][k]*dx_dxi[j][k]
        for k in range(3):
            b[i] += dx_dxi[i][k]*direction[k]
    # 2x2 matrix inverse
    deta = a[0][0]*a[1][1] - a[0][1]*a[1][0]
    if deta > 0.0:
        inva = [ [ a[1][1]/deta, -a[0][1]/deta ], [ -a[1][0]/deta, a[0][0]/deta ] ]
        delta_xi1 = inva[0][0]*b[0] + inva[0][1]*b[1]
        delta_xi2 = inva[1][0]*b[0] + inva[1][1]*b[1]
    else:
        # at pole: assume direction is inline with d1 or d2 and other is zero
        delta_xi2 = vector.dotproduct(d2, direction)
        if math.fabs(delta_xi2) > 0.0:
            delta_xi1 = 0.0
            delta_xi2 = (1.0 if (delta_xi2 > 0.0) else -1.0)*vector.magnitude(direction)/vector.magnitude(d2)
        else:
            delta_xi1 = vector.dotproduct(d1, direction)
            if math.fabs(delta_xi1) > 0.0:
                delta_xi1 = (1.0 if (delta_xi1 > 0.0) else -1.0)*vector.magnitude(direction)/vector.magnitude(d1)
                delta_xi2 = 0.0
    #delx = [ (delta_xi1*d1[c] + delta_xi2*d2[c]) for c in range(3) ]
    #print('delx', delx, 'dir', direction, 'diff', vector.magnitude([ (delx[c] - direction[c]) for c in range(3) ]))
    return delta_xi1, delta_xi2

def calculate_surface_axes(d1, d2, direction):
    '''
    :return: Vectors ax1, ax2, ax3: ax1 in-plane in 3-D vector direction,
    ax2 in-plane normal to a and ax3 normal to the surface plane.
    Vectors all have unit magnitude.
    '''
    delta_xi1, delta_xi2 = calculate_surface_delta_xi(d1, d2, direction)
    ax1 = vector.normalise([ delta_xi1*d1[c] + delta_xi2*d2[c] for c in range(3) ])
    ax3 = vector.crossproduct3(d1, d2)
    mag3 = vector.magnitude(ax3)
    if mag3 > 0.0:
        ax3 = [ s/mag3 for s in ax3 ]
        ax2 = vector.normalise(vector.crossproduct3(ax3, ax1))
    else:
        ax3 = [ 0.0, 0.0, 0.0 ]
        ax2 = [ 0.0, 0.0, 0.0 ]
    return ax1, ax2, ax3

def increment_xi_on_square(xi1, xi2, dxi1, dxi2):
    '''
    Increment xi1, xi2 by dxi1, dxi2 limited to square element bounds on [0,1].
    Works out face crossed first and limits that xi to 0.0 or 1.0 and other xi
    in proportion.
    :return: New xi1, xi2, proportion, face number 1-4 or None if within boundary.
    Proportion is 1.0 if increment is within element, otherwise equal to proportion
    of dxi1, dxi2 incremented to hit the boundary. If so limited to the boundary,
    the face number is returned as an integer:
    1 is xi1==0.0, 2 is xi1==1.0, 3 is xi2==0.0, 4 is xi2==1.0
    '''
    onxi1 = nxi1 = xi1 + dxi1
    onxi2 = nxi2 = xi2 + dxi2
    proportion = 1.0
    faceNumber = None
    if (nxi1 < 0.0) or (nxi1 > 1.0) or (nxi2 < 0.0) or (nxi2 > 1.0):
        # come back in direction of dxi1, dxi2 to first boundary
        if (onxi1 < 0.0) and (dxi1 < 0.0):
            thisProportion = -xi1/dxi1
            if thisProportion < proportion:
                proportion = thisProportion
                faceNumber = 1
                nxi1 = 0.0
                nxi2 = xi2 + proportion*dxi2
        elif (onxi1 > 1.0) and (dxi1 > 0.0):
            thisProportion = (1.0 - xi1)/dxi1
            if thisProportion < proportion:
                proportion = thisProportion
                faceNumber = 2
                nxi1 = 1.0
                nxi2 = xi2 + proportion*dxi2
        if (onxi2 < 0.0) and (dxi2 < 0.0):
            thisProportion = -xi2/dxi2
            if thisProportion < proportion:
                proportion = thisProportion
                faceNumber = 3
                nxi1 = xi1 + proportion*dxi1
                nxi2 = 0.0
        elif (onxi2 > 1.0) and (dxi2 > 0.0):
            thisProportion = (1.0 - xi2)/dxi2
            if thisProportion < proportion:
                proportion = thisProportion
                faceNumber = 4
                nxi1 = xi1 + proportion*dxi1
                nxi2 = 1.0
        #print('  increment to face', faceNumber, 'xi (' + str(xi1) + ',' + str(xi2) + ')', 'nxi (' + str(nxi1) + ',' + str(nxi2) + ')')
    return nxi1, nxi2, proportion, faceNumber
