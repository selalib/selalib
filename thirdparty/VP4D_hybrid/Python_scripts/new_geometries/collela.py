# -*- coding: UTF-8 -*-
"""
This scripts generates a Spline mapping that approachs the Collela
transformation
To run this script:
    python collela.py <nx> <ny> <px> <py>
"""

import numpy as np
from numpy import cos, sin, pi
import math
import matplotlib.pyplot as plt
import sys, os
from caid.cad_geometry import square

#-----------------------------------
try:
    nx = int(sys.argv[1])
except:
    nx = 3

try:
    ny = int(sys.argv[2])
except:
    ny = 3

try:
    px = int(sys.argv[3])
except:
    px = 3

try:
    py = int(sys.argv[4])
except:
    py = 3
#-----------------------------------

# ...
geo = square(n=[nx, ny], p=[px,py])
# ...

# ...
nrb = geo[0]

u = nrb.greville(axis=0)
v = nrb.greville(axis=1)

U,V = np.meshgrid(u,v)

u = U.transpose()
v = V.transpose()
# ...

# ...
eps = 0.1
k1 = 1. ; k2 = 1.
X = 2*(u + eps * sin(2*pi*k1*u) * sin(2*pi*k2*v)) -1.
Y = 2*(v + eps * sin(2*pi*k1*u) * sin(2*pi*k2*v)) - 1.
# ...

# ...

ctrl_pts = np.zeros_like(geo[0].points)
ctrl_pts[:,:,0] = X
ctrl_pts[:,:,1] = Y
geo[0].set_points(ctrl_pts)
# ...

geo.plotMesh()
plt.show()

# ... export to txt format
geo.save("collela.txt")
# ...

