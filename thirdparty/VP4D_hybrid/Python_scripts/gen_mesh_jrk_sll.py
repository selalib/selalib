#! /usr/bin/python

from caid.cad_geometry import square
from caid.io import NML
from numpy import pi
import os
import sys

# -----------------------------------------------
try:
    nx = int(sys.argv[1])
    ny = int(sys.argv[2])
except:
    print ("you must run this script with python gen_mesh_caid_custom_rescale.py n1 n2.")
    print ("n1 and n2 denote the numer of elements in the x and y directions")
    raise()
# -----------------------------------------------


geometry_path = "geometry/"
os.system("mkdir -p " + geometry_path)

L = [4. * pi , 4. * pi]
geo = square(n=[nx-1,ny-1], p=[3,3])

geo._internal_faces = []
geo._external_faces = []
geo._connectivity   = []

dict_con = {}
dict_con['original'] = [0,0]; dict_con['clone'] = [0,2]
geo._connectivity.append(dict_con)
dict_con = {}
dict_con['original'] = [0,1]; dict_con['clone'] = [0,3]
geo._connectivity.append(dict_con)

for axis in range(0, 2):
    geo.scale(L[axis], axis=axis)

# ... export JOREK format
list_NodeData, list_ElementData= geo.to_bezier_jorek(0, "jorek")

os.system("mv jorek_elements.txt " + geometry_path)
os.system("mv jorek_nodes.txt " + geometry_path)
# ...

# ... export selalib format
selalib_geometry_path = 'sll_geometry.nml'

nml_io = NML()
nml_io.write(selalib_geometry_path, geo)

os.system("mv sll_geometry_*.nml " + geometry_path)
# ...
