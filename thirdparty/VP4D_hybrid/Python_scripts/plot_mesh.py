import matplotlib.pyplot as mpp
import numpy as np
import os
import VP4Dhyb_GYSutils
GYSut = VP4Dhyb_GYSutils.impall('VP4Dhyb_GYSutils')
VP4Dd = GYSut.impall('VP4Dhyb_diag')

mesh_directory = raw_input("Directory for mesh [default = current] ?:");
if (mesh_directory==''):
    mesh_directory = str(os.getcwd())

ff = GYSut.loadHDF5(mesh_directory+'/mesh.h5')
VP4Dd.mesh_plot(ff)
