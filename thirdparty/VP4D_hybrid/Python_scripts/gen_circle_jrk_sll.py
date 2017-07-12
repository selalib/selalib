#! /usr/bin/python

from caid.cad_geometry import circle, annulus
from caid.io import NML
from numpy import pi
import os
import sys

# -----------------------------------------------
try:
    nx = int(sys.argv[1])
    ny = int(sys.argv[2])
    print ("Creating annulus geometry with nx, ny =", nx, ny)
except:
    print ("you must run this script with python gen_circle_jrk_sll.py n1 n2.")
    print ("n1 and n2 denote the numer of elements in the x and y directions")
    raise()
# -----------------------------------------------


#geo = circle(n=[nx-1,ny-1], center=[0.,0.], radius=4.0*pi)
geo = annulus(n=[nx,ny], center=[0.,0.], rmin=0.1, rmax=12.0, p=[3, 3])
# geo.save("jorek")


# ... export JOREK format
list_NodeData, list_ElementData= geo.to_bezier_jorek(0, "jorek")
# ...

# ... export selalib format
selalib_geometry_path = 'sll_geometry.nml'

nml_io = NML()
nml_io.write(selalib_geometry_path, geo)
# ...
