'''
Utility functions for vectors.
'''

import math

def crossproduct3(a, b):
    '''
    :return: vector 3-D cross product of a and b
    '''
    return [ a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0] ]

def dotproduct(a, b):
    '''
    :return: vector dot (inner) product of a and b
    '''
    return sum(a[i]*b[i] for i in range(len(a)))

def magnitude(v):
    '''
    return: scalar magnitude of vector v
    '''
    return math.sqrt(sum(c*c for c in v))

def normalise(v):
    '''
    :return: vector v normalised to unit length
    '''
    mag = math.sqrt(sum(c*c for c in v))
    return [ c/mag for c in v ]

def setMagnitude(v, mag):
    '''
    return: Vector v with magnitude set to mag.
    '''
    scale = mag/math.sqrt(sum(c*c for c in v))
    return [ c*scale for c in v ]

def addVectors(v1,v2,s1=1,s2=1):
    '''
    returns s1*v1+s2*v2 where s1 and s2 are scalars.
    :return: Vector s1*v1+s2*v2
    '''
    return [s1 * v1[c] + s2 * v2[c] for c in range(len(v1))]
