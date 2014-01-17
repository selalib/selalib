import utils_func as ut
import matplotlib as mp
import matplotlib.pyplot as mpp

filename = '/Users/back/selalib/prototype/build/DK4d_diag_d00002.h5'
print "HDF5 file read = " + filename
f_diag = ut.loadHDF5(filename)

fig = mpp.figure(figsize=(18,9))

ax1 = fig.add_subplot(1,5,5)
p1  = ax1.pcolormesh(f_diag.xgrid_2d,f_diag.ygrid_2d,f_diag.phi2d_xy)
mpp.xlabel('x')
mpp.ylabel('y')
mpp.title('Phi2D_xy')
fig.colorbar(p1)

ax2 = fig.add_subplot(1,5,3)
p2  = ax2.pcolormesh(f_diag.xgrid_2d,f_diag.ygrid_2d,f_diag.rho2d_xy)
mpp.xlabel('x')
mpp.ylabel('y')
mpp.title('rho2D_xy')
fig.colorbar(p2)

ax3 = fig.add_subplot(1,5,4)
p3  = ax3.pcolormesh(f_diag.xgrid_2d,f_diag.ygrid_2d,f_diag.f2d_xy)
mpp.xlabel('x')
mpp.ylabel('y')
mpp.title('f2D_xy')
fig.colorbar(p3)

#ax4 = fig.add_subplot(2,6,4)
#p4  = ax4.pcolormesh(f_diag.f2d_zvpar)
#mpp.xlabel('z')
#mpp.ylabel('vpar')#
#mpp.title('f2D_zvpar')
#fig.colorbar(p4)

ax5 = fig.add_subplot(1,5,1)
p5  = ax5.pcolormesh(f_diag.xgrid_2d,f_diag.ygrid_2d,f_diag.E2d_eta1_xy)
mpp.xlabel('x')
mpp.ylabel('y')
mpp.title('E2d_x1_xy')
fig.colorbar(p5)

ax6 = fig.add_subplot(1,5,2)
p6  = ax6.pcolormesh(f_diag.xgrid_2d,f_diag.ygrid_2d,f_diag.E2d_eta2_xy)
mpp.xlabel('x')
mpp.ylabel('y')
mpp.title('E2d_x2_xy')
fig.colorbar(p6)


#ax7 = fig.add_subplot(2,6,6)
#mpp.xlabel('x')
#mpp.ylabel('y')
#mpp.title('nrj_kin')
#mpp.plot(f_diag.nrj_kin)

#ax8 = fig.add_subplot(2,6,7)
#mpp.xlabel('x')
#mpp.ylabel('y')
#mpp.title('nrj_pot')
#mpp.plot(f_diag.nrj_pot)

#ax9 = fig.add_subplot(2,6,8)
#mpp.xlabel('x')
#mpp.ylabel('y')
#mpp.title('nrj_tot')
#mpp.plot(f_diag.nrj_tot)



#ax10 = fig.add_subplot(2,6,9)
#mpp.xlabel('x')
#mpp.ylabel('y')

#mpp.hold(True)
#mpp.plot(f_diag.r_grid,f_diag.n0_r,'r+')
#mpp.plot(f_diag.r_grid,f_diag.n0_r_func,'b')
#mpp.title('profil n0_r')
#mpp.hold(False)

#ax11 = fig.add_subplot(2,6,10)
#mpp.xlabel('x')
#mpp.ylabel('y')
#mpp.hold(True)
#mpp.plot(f_diag.r_grid,f_diag.Ti_r,'r+')
#mpp.plot(f_diag.r_grid,f_diag.Ti_r_func,'b')
#mpp.title('profil Ti_r')
#mpp.hold(False)

#ax12 = fig.add_subplot(2,6,11)
#mpp.xlabel('x')
#mpp.ylabel('y')

#mpp.hold(True)
#mpp.plot(f_diag.r_grid,f_diag.Te_r,'r+')
#mpp.plot(f_diag.r_grid,f_diag.Te_r_func,'b')
#mpp.title('profil Te_r')
#mpp.hold(False)

#ax10 = fig.add_subplot(2,6,9)
#mpp.xlabel('x')
#mpp.ylabel('y')
#mpp.title('relative_error_masse')
#mpp.plot(f_diag.relative_error_masse)

#ax11 = fig.add_subplot(2,6,10)
#mpp.xlabel('x')
#mpp.ylabel('y')
#mpp.title('L1_norm')
#mpp.plot(f_diag.L1_norm)

#ax12 = fig.add_subplot(2,6,11)
#mpp.xlabel('x')
#mpp.ylabel('y')
#mpp.title('L2_norm')
#mpp.plot(f_diag.L2_norm)

#ax13 = fig.add_subplot(2,6,12)
#mpp.xlabel('x')
#mpp.ylabel('y')
#mpp.title('Linf_norm')
#mpp.plot(f_diag.Linf_norm)

#fig1=mpp.figure()
#mpp.plot(f_diag.xgrid_2d,f_diag.ygrid_2d,'+')
#fig1.show()

fig.show()


#fig2=mpp.figure()
#mpp.plot(f_diag.xgrid_2d,f_diag.ygrid_2d,f_diag.rho2d_xy,'+')
#fig2.show()

#fig3=mpp.figure()
#mpp.plot(f_diag.xgrid_2d,f_diag.ygrid_2d,f_diag.f2d_zvpar,'-')
#fig3.show()

