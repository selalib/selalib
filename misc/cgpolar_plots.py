import glob
import h5py
import numpy as np
import pylab as pl



#mesh
f = h5py.File("polar_mesh-x1.h5","r")
x1 = f['x1'].value
f.close()
f = h5py.File("polar_mesh-x2.h5","r")
x2 = f['x2'].value
f.close()

nfile = len(glob.glob1('.',"f*-values.h5"))


