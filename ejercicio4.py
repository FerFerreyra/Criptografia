p = 101
a = 13
b = 37

# Create a simple Point class to represent the affine points.
from collections import namedtuple
Point = namedtuple("Point", "x y")

# The point at infinity (origin for the group law).
O = 'Origin'

def valid(P):
    """
    Determine whether we have a valid representation of a point
    on our curve.  We assume that the x and y coordinates
    are always reduced modulo p, so that we can compare
    two points for equality with a simple ==.
    """
    if P == O:
        return True
    else:
        return (
            (P.y**2 - (P.x**3 + a*P.x + b)) % p == 0 and
            0 <= P.x < p and 0 <= P.y < p)

def inv_mod_p(x):
    """
    Compute an inverse for x modulo p, assuming that x
    is not divisible by p.
    """
    if x % p == 0:
        raise ZeroDivisionError("Impossible inverse")
    return pow(x, p-2, p)

def ec_inv(P):
    """
    Inverse of the point P on the elliptic curve y^2 = x^3 + ax + b.
    """
    if P == O:
        return P
    return Point(P.x, (-P.y)%p)

def ec_add(P, Q):
    """
    Sum of the points P and Q on the elliptic curve y^2 = x^3 + ax + b.
    """
    if not (valid(P) and valid(Q)):
        raise ValueError("Invalid inputs")

    # Deal with the special cases where either P, Q, or P + Q is
    # the origin.
    if P == O:
        result = Q
    elif Q == O:
        result = P
    elif Q == ec_inv(P):
        result = O
    else:
        # Cases not involving the origin.
        if P == Q:
            dydx = (3 * P.x**2 + a) * inv_mod_p(2 * P.y)
        else:
            dydx = (Q.y - P.y) * inv_mod_p(Q.x - P.x)
        x = (dydx**2 - P.x - Q.x) % p
        y = (dydx * (P.x - x) - P.y) % p
        result = Point(x, y)

    # The above computations *should* have given us another point
    # on the curve.
    assert valid(result)
    return result

P = Point(36,87)
Q = P 
print(ec_add(P,P))

def ec_mul(P, n):
    if (n == 0) :
        return O
    elif (n == 1) :
        return P
    else : 
        Q = O
        for i in range (n):
            Q = ec_add(Q,P)
        return Q

print(ec_mul(P,108))

mensaje = [24,0,16,20,8,4,17,14,21,0,2,0,2,8,14,13,4,18]

def num2point(P,m):
    points = []
    for i in range(len(m)):
        points = points + [ec_mul(P,m[i])]
    return points

mensaje = num2point(P,mensaje)

x1=53
x2=87
y1= ec_mul(P,x1)
y2= ec_mul(P,x2)
k=37

print(y1)
print(y2)

def gamal(P,mensaje,key,k):
    cifrado = []
    for m in mensaje: 
        cifrado = cifrado + [ec_add(m,ec_mul(key,k))]
    return (ec_mul(P,k), cifrado)

print(gamal(P,mensaje,y1,k))
