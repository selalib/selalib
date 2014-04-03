import utils_func as ut
import matplotlib as mp
import matplotlib.pyplot as mpp

filename = '/Users/back/selalib/prototype/build/vp4D_diag_d00004.h5'
f_diag = ut.loadHDF5(filename)
print "HDF5 file read = " + filename


fig = mpp.figure(figsize=(18,9))

ax3 = fig.add_subplot(1,2,1)
p3  = ax3.pcolormesh(f_diag.f2d_xy)
mpp.xlabel('x')
mpp.ylabel('y')
mpp.title('f2D_xy')
fig.colorbar(p3)

ax2 = fig.add_subplot(1,2,2)
p2  = ax2.pcolormesh(f_diag.f2d_v1v2)
mpp.xlabel('vx')
mpp.ylabel('vy')
mpp.title('f2D_v1v2')
fig.colorbar(p2)

fig.show()
