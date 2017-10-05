# coding: utf-8
import matplotlib.pylab as mpp
from caid.cad_geometry import square
import numpy as np
from caid.io import NML

n = 60
p = 3 
alpha_coeff = 0.95  #0.95

# Create a first coarse mesh
geo = square( p=[p,p] )
patch0 = geo[0] 
pts = patch0.points

# Modify the position of the control points
regular = 1./3.
alpha = alpha_coeff*regular
pts[1,1,:2] =  [regular-alpha,regular-alpha]       #regular#1./3.
pts[1,2,:2] =  [regular-alpha,2.*regular+alpha]    #regular#1./3.
pts[2,1,:2] =  [2.*regular+alpha,regular-alpha]    #regular#1./3.
pts[2,2,:2] =  [2.*regular+alpha,2.*regular+alpha] #regular#1./3.
print "1 = ", pts[1,1]
print "2 = ", pts[1,2]
print "3 = ", pts[2,1]
print "4 = ", pts[2,2]
patch0.set_points(pts)

# Create the finer mesh
t = np.linspace(0.,1.,n+2)[1:-1] # remove 0 and 1, only internal knots will be inserted
geo.refine(list_t=[t,t])

# Rescale the square to 4pi
geo.scale(4.*np.pi)

#---> Plot mesh
mpp.figure(figsize=(6,6))
geo.plotMesh(MeshResolution=5)
#-----> Save the geometry
s_alpha_coeff = str(alpha_coeff).replace('.','')
s_alpha_coeff = s_alpha_coeff.replace('-','m')
caid_acronym = "square_4pi_modif_"+s_alpha_coeff+"_n"+str(n)+"n"+str(n)
geo.save(caid_acronym+".txt")
nml_io = NML()
nml_io.write(caid_acronym+".nml",geo)
mpp.savefig("MESH_"+caid_acronym+".png")
mpp.show()

#---> Plot jacobian
mpp.figure(figsize=(6,6))
geo.plotJacobians(MeshResolution=50)
mpp.savefig("JACOB_"+caid_acronym+".png")
mpp.colorbar()
mpp.show()




