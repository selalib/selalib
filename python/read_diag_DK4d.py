import utils_func as ut
import matplotlib as mp
import matplotlib.pyplot as mpp

filename = 'DK4d_diag_d00000.h5'
f_diag = ut.loadHDF5(filename)
print "HDF5 file read = " + filename

fig = mpp.figure(figsize=(12,12))

ax1 = fig.add_subplot(2,2,1)
p1  = ax1.pcolormesh(f_diag.phi2d_xy)
mpp.xlabel('y')
mpp.ylabel('x')
mpp.title('Phi2D_xy')
fig.colorbar(p1)

ax2 = fig.add_subplot(2,2,2)
p2  = ax2.pcolormesh(f_diag.rho2d_xy)
mpp.xlabel('y')
mpp.ylabel('x')
mpp.title('rho2D_xy')
fig.colorbar(p2)

ax3 = fig.add_subplot(2,2,3)
p3  = ax3.pcolormesh(f_diag.f2d_xy)
mpp.xlabel('y')
mpp.ylabel('x')
mpp.title('f2D_xy')
fig.colorbar(p3)

ax4 = fig.add_subplot(2,2,4)
p4  = ax4.pcolormesh(f_diag.f2d_zvpar)
mpp.xlabel('vpar')
mpp.ylabel('z')
mpp.title('f2D_zvpar')
fig.colorbar(p4)

fig.show()


