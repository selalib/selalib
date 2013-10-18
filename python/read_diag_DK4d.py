import utils_func as ut
import matplotlib as mp
import matplotlib.pyplot as mpp

f_diag = ut.loadHDF5('/Users/back/selalib/prototype/build/DK4d_diag_d00000.h5')

fig = mpp.figure(figsize=(10,8))

ax1 = fig.add_subplot(2,2,1)
p1  = ax1.pcolormesh(f_diag.phi2d_x1x2)
mpp.xlabel('x2')
mpp.ylabel('x1')
mpp.title('Phi2D_x1x2')
fig.colorbar(p1)

ax2 = fig.add_subplot(2,2,2)
p2  = ax2.pcolormesh(f_diag.rho2d_x1x2)
mpp.xlabel('x2')
mpp.ylabel('x1')
mpp.title('rho2D_x1x2')
fig.colorbar(p2)

ax3 = fig.add_subplot(2,2,3)
p3  = ax3.pcolormesh(f_diag.f2d_x1x2)
mpp.xlabel('x2')
mpp.ylabel('x1')
mpp.title('f2D_x1x2')
fig.colorbar(p3)

ax4 = fig.add_subplot(2,2,4)
p4  = ax4.pcolormesh(f_diag.f2d_x3x4)
mpp.xlabel('x4')
mpp.ylabel('x3')
mpp.title('f2D_x3x4')
fig.colorbar(p4)

fig.show()


