import utils_func as ut
import matplotlib as mp
import matplotlib.pyplot as mpp

f_diag = ut.loadHDF5('/Users/back/selalib/prototype/build/DK4d_diag_d00010.h5')

fig = mpp.figure(figsize=(14,8))

ax1 = fig.add_subplot(2,3,1)
p1  = ax1.pcolormesh(f_diag.phi2d_xy)
mpp.xlabel('x')
mpp.ylabel('y')
mpp.title('Phi2D_xy')
fig.colorbar(p1)

ax2 = fig.add_subplot(2,3,2)
p2  = ax2.pcolormesh(f_diag.rho2d_xy)
mpp.xlabel('x')
mpp.ylabel('y')
mpp.title('rho2D_xy')
fig.colorbar(p2)

ax3 = fig.add_subplot(2,3,3)
p3  = ax3.pcolormesh(f_diag.f2d_xy)
mpp.xlabel('x')
mpp.ylabel('y')
mpp.title('f2D_xy')
fig.colorbar(p3)

ax4 = fig.add_subplot(2,3,4)
p4  = ax4.pcolormesh(f_diag.f2d_zvpar)
mpp.xlabel('z')
mpp.ylabel('vpar')#
mpp.title('f2D_zvpar')
fig.colorbar(p4)

ax5 = fig.add_subplot(2,3,5)
p5  = ax5.pcolormesh(f_diag.E2d_eta1_xy)
mpp.xlabel('x')
mpp.ylabel('y')
mpp.title('E2d_eta1_xy')
fig.colorbar(p5)

ax6 = fig.add_subplot(2,3,6)
p6  = ax6.pcolormesh(f_diag.E2d_eta2_xy)
mpp.xlabel('x')
mpp.ylabel('y')
mpp.title('E2d_eta2_xy')
fig.colorbar(p6)


fig.show()


