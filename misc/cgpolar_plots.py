import glob
import h5py
import matplotlib.pyplot as pl

h5file = h5py.File("polar_mesh-x1.h5","r")
x1 = h5file['x1'].value
h5file.close()
h5file = h5py.File("polar_mesh-x2.h5","r")
x2 = h5file['x2'].value
h5file.close()

nfile = len(glob.glob1('.',"f*-values.h5"))-1
h5filename = "f%04d-values.h5" % (nfile)
print h5filename
h5file = h5py.File(h5filename,"r")
f = h5file['values'].value
h5file.close()

surf = pl.contourf(x1,x2,f)
cb = pl.colorbar(surf, shrink=0.5, aspect=10)
pl.show()
