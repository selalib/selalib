import math

#--------------------------------------------
#
# Class to define a point with cartesian coordinates
#
#--------------------------------------------
class Point(object):

    def __init__(self, x, y):
        self.X = x
        self.Y = y

    def translation(self, dx, dy):
        self.X = self.X + dx
        self.Y = self.Y + dy

    def distance(self, other):
        # Distance to origin (0,0)
        dx = self.X - other.X
        dy = self.Y - other.Y
        return math.hypot(dx,dy)

    def __add__(self,other):
        x = self.X + other.X
        y = self.Y + other.Y
        return Point(x,y)

    def __mul__(self, val):
        return Point(self.X * val, self.Y * val)

    def __eq__(self, other):
        return ((self.X==other.X)and(self.Y==other.Y))

    def __neq__(self, other):
        return ((self.X<>other.X)or(self.Y<>other.Y))

    def __str__(self):
        return "(%s,%s)"%(self.X, self.Y)


class Vector(Point):
    def __init__(self, p1, p2):
        self.X = p2.X - p1.X
        self.Y = p2.Y - p1.Y

    def __str__(self):
        return "(%s,%s)"%(self.X, self.Y)

    # def __add__(self, other):
    #     res = Vector(self.X + other.X, self.Y + other.Y)
    #     return res

    # def __mul__(self, other):
    #     res = Vector(self.X * other.X, self.Y * other.Y)
    #     return res

    # def __mul__(self, val):
    #     res = Vector(self.X * val, self.Y * val)
    #     return res

    def norm(self, other):
        # Same as distance for a point
        dx = self.X - other.X
        dy = self.Y - other.Y
        return math.hypot(dx,dy)


#**********************************************************
# Defines and returns coordinates of the Fekete
# points of degree 3, plus associated weights
# Reference:
#      Mark Taylor, Beth Wingate, Rachel Vincent,
#      An Algorithm for Computing Fekete Points in the Triangle,
#      SIAM Journal on Numerical Analysis,
#      Volume 38, Number 5, 2000, pages 1707-1720.
#**********************************************************
def fekete3(p1, p2, p3) :
    # Coordinates of base orbits points as defined in Ref
    orbit1  = Point(1./3., 1./3.)
    orbit2  = Point(0.0  , 0.0  )
    orbit3  = Point(0.0  , 0.2763932023)
    orbits  = [orbit1, orbit2, orbit3]
    weight1 = 0.45
    weight2 = 1./60.
    weight3 = 1./12.
    weights = [weight1, weight2, weight3]

    # List containing points
    fekPts  = []
    fekWei  = []

    # Defining Fekete points parting from p1
    v1 = Vector(p1, p2)
    v2 = Vector(p1, p3)
    for i in range(3) :
        orb = orbits[i]
        new_fek  = (v1*orb.X + v2*orb.Y) + p1
        new_fek2 = (v1*orb.Y + v2*orb.X) + p1
        if (fekPts.count(new_fek) == 0) :
            fekPts.append(new_fek)
            fekWei.append(weights[i])
        if (fekPts.count(new_fek2) == 0) :
            fekPts.append(new_fek2)
            fekWei.append(weights[i])

    # Defining Fekete points parting from p2
    v1 = Vector(p2, p1)
    v2 = Vector(p2, p3)
    for i in range(1,3) :
        orb = orbits[i]
        new_fek  = (v1*orb.X + v2*orb.Y) + p2
        new_fek2 = (v1*orb.Y + v2*orb.X) + p2
        if (fekPts.count(new_fek) == 0) :
            fekPts.append(new_fek)
            fekWei.append(weights[i])
        if (fekPts.count(new_fek2) == 0) :
            fekPts.append(new_fek2)
            fekWei.append(weights[i])

    # Defining Fekete points parting from p3
    v1 = Vector(p3, p2)
    v2 = Vector(p3, p1)
    for i in range(1,3) :
        orb = orbits[i]
        new_fek  = (v1*orb.X + v2*orb.Y) + p3
        new_fek2 = (v1*orb.Y + v2*orb.X) + p3
        if (fekPts.count(new_fek) == 0) :
            fekPts.append(new_fek)
            fekWei.append(weights[i])
        if (fekPts.count(new_fek2) == 0) :
            fekPts.append(new_fek2)
            fekWei.append(weights[i])

    return [fekPts, fekWei]


def main():
    #Defining the vertices of the triangle
    p1 = Point(0,0)
    p2 = Point(1,0)
    p3 = Point(0,1)

    [fekPts, fekWei] = fekete3(p1, p2, p3)

    for i in range(len(fekPts)) :
        print '{0:2d} {1:6f} {2:6f} {3:6f}'.format(i+1, fekPts[i].X, fekPts[i].Y, fekWei[i])

main()
