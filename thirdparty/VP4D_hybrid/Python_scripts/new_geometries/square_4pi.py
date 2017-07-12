# coding: utf-8
import matplotlib.pylab as mpp
from caid.cad_geometry import square

mpp.ion()
geo = square(n=[63,63], p=[3,3])
from numpy import pi
geo.scale(4.*pi)
geo[0].points
geo.plotMesh()
geo.save("square_4pi_n63n63p3p3.txt")
mpp.ioff()
