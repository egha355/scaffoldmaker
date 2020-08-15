'''
Utilities for getting a mirror image of a scaffold.
'''

from __future__ import division
from scaffoldmaker.utils import vector

def getPointDistanceFromPlane(x, p):
    """
    Get distance from point x to plane p. Returns negative distance, if point x is on the negative side of the plane.
    :param x: Point coordinates.
    :param p: Plane ax+by+cd=d defined as a list p=[a,b,c,d].
    :return: Distance.
    """
    n = getPlaneNormalVector(p)
    ma = vector.magnitude(p[:2])
    return vector.dotproduct(x, n) - p[3] / ma

def getPlaneNormalVector(p):
    """
    Get plane normal vector
    :param p:  Plane ax+by+cd=d defined as a list p=[a,b,c,d].
    :return: Unit normal vector.
    """
    return vector.normalise([p[0], p[1], p[2]])

def pointProjectionToPlane(x, p):
    """
    Find the projection of point x onto plane p
    :param x:  Point coordinates.
    :param p: Plane ax+by+cd=d defined as a list p=[a,b,c,d].
    :return: x0, the point projection onto the plane
    """
    D = getPointDistanceFromPlane(x, p)
    n = getPlaneNormalVector(p)
    return [x[c] - D*n[c] for c in range(3)]

def mirrorImageOfPoint(x, p):
    """
    Find the mirror image of point x on plane p
    :param x:  Point coordinates.
    :param p: Plane ax+by+cd=d defined as a list p=[a,b,c,d].
    :return: xm, mirror image of point x.
    """
    D = getPointDistanceFromPlane(x, p)
    n = getPlaneNormalVector(p)
    return [x[c] - 2*D*n[c] for c in range(3)]

def mirrorVector(v, p, keepSign = True):
    """
    Find the mirror image of vector v on plane p
    :param v: Vector.
    :param p: Plane.
    :param keepSign: If True, negative of the vector is returned.
    :return: image vector.
    """
    s = 1 if keepSign else -1
    n = getPlaneNormalVector(p)
    return [s*v[c] - s*2 * vector.dotproduct(v, n) * n[c] for c in range(3)]