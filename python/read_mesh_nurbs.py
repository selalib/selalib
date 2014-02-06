import utils_func as ut
import matplotlib as mp
import matplotlib.pyplot as mpp

filename_x1 = '/Users/back/selalib/prototype/build/mhdeq_patch0-x1.h5'
x1 = ut.loadHDF5(filename_x1)
print "HDF5 file read = " + filename_x1

filename_x2 = '/Users/back/selalib/prototype/build/mhdeq_patch0-x2.h5'
x2 = ut.loadHDF5(filename_x2)
print "HDF5 file read = " + filename_x2


fig=mpp.figure()

mpp.plot(x1.x1,x2.x2,'+')
mpp.plot(x2.x2,x1.x1,'+')
fig.show()
