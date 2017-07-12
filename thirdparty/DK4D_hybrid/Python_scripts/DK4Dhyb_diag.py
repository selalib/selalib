#-------------------------------------------------
#  file : DK4Dhyb_diag.py
#  date : 2014-08-25
#
#   Diagnostics for DK4D code
#-------------------------------------------------
import GYSutils
import matplotlib.pyplot as mpp
import numpy as np
GYSut      = GYSutils.impall('GYSutils')

#-----------------------------------------------
# Diagnostic choices
#-----------------------------------------------
def DK4Ddiag(resu):
    """ Diagnostic choices"""

    choice = '1';
    while choice != '' :
        print ' '
        print ' **** EQUILIBRIUM ****'
        print '  *   1 : (x,y) density and temperature profiles'
        print '  *   2 : equilibrium distribution function'
        print ' '
        print ' **** PHI2D DIAGNOSTICS ****'
        print '  *   11 : Phi(x,y)'
        print '  *   12 : int Phi(x,y,z)**2 jac dx dy dz'
        print ' '

        print ' **** F2D DIAGNOSTICS ****'
        print '  *   21 : delta f2D(x,y)'
        print '  *   22 : delta f2D(z,vpar)'
        print ' '

        print ' **** CONSERVATION LAWS ****'
        print '  *   31 : L1-Norm, L2-Norm and Linf-Norm'
        print '  *   32 : Energy'
        print ' '

        choice = raw_input(' choice ? : ');
        choice = GYSut.str2num(choice)

        if (choice != ''):
            DK4Ddiag_choice(choice,resu)

    #end while

    print 'Bye bye!'

#end def DK4Ddiag


#-----------------------------------------------------------------
#  Execution of DK4D diagnostics according to the user choice
#-----------------------------------------------------------------
def DK4Ddiag_choice(ichoice,resu):
    """
    Execution of GYSELA diagnostics according to the user choice
    """

    # Phi2D diagnostics
    if (ichoice == 1):
        equilibrium_profile_plot(resu)
    if (ichoice == 2):
        equilibrium_fdistribu_plot(resu)
    elif (ichoice == 11):
        Phi2D_xy_plot(resu)
    elif (ichoice == 12):
        Phisquare_plot(resu)
    elif (ichoice == 21):
        deltaf2D_xy_plot(resu)
    elif (ichoice == 22):
        deltaf2D_eta3vpar_plot(resu)
    elif (ichoice == 31):
        L1_L2_Linf_norms_plot(resu)
    elif (ichoice == 32):
        energies_plot(resu)

#end def DK4Ddiag_choice


#------------------------------------------------------
# Plot all the quantities related to the equilibrium
# --------
# choice 1
#------------------------------------------------------
def equilibrium_profile_plot(resu):

    """ Equilibrium profile"""

    fig = mpp.figure(figsize=(6,6))

    ax1 = fig.add_subplot(3,2,1)
    mpp.plot(resu.rgrid_2d[0][:],resu.Ti_r,'-+',linewidth=2) 
    mpp.xlabel('r')
    mpp.ylabel('T_i')
    mpp.grid(True)
    
    ax2 = fig.add_subplot(3,2,2)
    mpp.plot(resu.rgrid_2d[0][:],resu.Te_r,'-+',linewidth=2)
    mpp.xlabel('r')
    mpp.ylabel('T_e')
    mpp.grid(True)
    
    ax3 = fig.add_subplot(3,2,3)
    mpp.plot(resu.rgrid_2d[0][:],resu.n0_r,'-+',linewidth=2)
    mpp.xlabel('r')
    mpp.ylabel('n_0')
    mpp.grid(True)
        
    Nr = np.size(resu.eta1_grid)-1
    ax4 = fig.add_subplot(3,2,4)
    ns0onTs0 = resu.n0_r[0:Nr]/resu.Ti_r[0:Nr]
    diff_ns0 = np.diff(resu.n0_r)
    diff_ns0[np.where(diff_ns0 == 0)] = 1.e-15
    eta = ns0onTs0*np.diff(resu.Ti_r)/diff_ns0
    mpp.plot(resu.rgrid_2d[0][0:Nr],eta,'-+',linewidth=2)
    mpp.xlabel('r')
    mpp.ylabel('\eta')
    mpp.grid(True)

    ax5 = fig.add_subplot(3,2,5)
    p5 = ax5.pcolor(resu.xgrid_2d,resu.ygrid_2d,resu.Ti_xy)
    fig.colorbar(p5)
    mpp.xlabel('x')
    mpp.ylabel('y')
    mpp.title('Ti(x,y)')
    mpp.grid(True)

    ax6 = fig.add_subplot(3,2,6)
    p6 = ax6.pcolor(resu.xgrid_2d,resu.ygrid_2d,resu.n0_xy)
    fig.colorbar(p6)
    mpp.xlabel('x')
    mpp.ylabel('y')
    mpp.title('n0(x,y)')
    mpp.grid(True)

#end def equilibrium_profile_plot


#------------------------------------------------------
# Plot the equilibrium distribution function
# --------
# choice 2
#------------------------------------------------------
def equilibrium_fdistribu_plot(resu):

    """ Equilibrium distribution function """

    fig = mpp.figure(figsize=(6,6))

    indx_x1 = resu.fdistribu_indx_diag[0][0]
    indx_x2 = resu.fdistribu_indx_diag[0][1]
    indx_x4 = resu.fdistribu_indx_diag[0][3]

    ax1 = fig.add_subplot(2,1,1)
    p1 = ax1.pcolor( 
        resu.xgrid_2d,
        resu.ygrid_2d,
        resu.feq_xyvpar[indx_x4,:,:] )
    fig.colorbar(p1)
    mpp.xlabel('x')
    mpp.ylabel('y')
    mpp.title('feq(x,y)')
    mpp.grid(True)

    ax2 = fig.add_subplot(2,1,2)
    p2 = ax2.plot(
            resu.vpar_grid,
            resu.feq_xyvpar[:,indx_x2,indx_x1],
            )
    mpp.xlabel('vpar')
    mpp.ylabel('feq(vpar)')
    mpp.grid(True)

#end def equilibrium_profile_plot


#------------------------------------------------------
# Phi(x,y)
# --------
# choice 11
#------------------------------------------------------
def Phi2D_xy_plot(resu):
    """ Phi(x,y)"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.diag_time,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = ax1.pcolormesh( resu.xgrid_2d[:,:], 
                              resu.ygrid_2d[:,:], 
                              resu.Phi_x1x2[itime_diag,:,:] ) 

        fig.colorbar(p1)
        fig.show()

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)
    #end while

#end def Phi2D_xy_plot



#------------------------------------------------------
# log int Phi(eta1,eta2,eta3)**2 jac deta1 deta2 deta3
# --------
# choice 12
#------------------------------------------------------
def Phisquare_plot(resu):
    """ Phi(x,y)"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.diag_time,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = mpp.semilogy( 
            resu.diag_time[0:itime_diag + 1], 
            resu.Phisquare[0:itime_diag + 1],'+-' )
        mpp.xlabel('time')
        mpp.ylabel('int Phi^2')
        mpp.grid(True)
        fig.show()

        #--> Computation of the growth rate
        print ' ##### parameters for the linear interpolation #####'
        [stime_init_interp,itime_init_interp] = GYSut.Ask_time(resu.diag_time, \
               ask_msg='Time for interpolation beginning',give_default=1)
        [stime_end_interp,itime_end_interp] = GYSut.Ask_time(resu.diag_time, \
               first_indx_time=itime_init_interp,               
               ask_msg='Time for interpolation end',give_default=2)
        time_tmp      = resu.diag_time[itime_init_interp:itime_end_interp+1].T[0]
        Phisquare_tmp = resu.Phisquare[itime_init_interp:itime_end_interp+1].T[0]
        poly_resu = np.polyfit( time_tmp, np.log(Phisquare_tmp), 1, full='True' )
        [p,s]           = poly_resu[0]
        residuals       = poly_resu[1]
        rank            = poly_resu[2] 
        singular_values = poly_resu[3]
        rcond           = poly_resu[4]
        print '=> growth_rate = ', p
        print '=> residuals   = ', residuals

        mpp.hold(True)
        str_legend = r'$\gamma = $' + str(p) + '\n residual = ' + str(residuals)
        mpp.semilogy( resu.diag_time[0:itime_end_interp+1], 
                      np.exp(s)*np.exp(p*resu.diag_time[0:itime_end_interp+1]),'r',
                      label=str_legend )
        mpp.hold(False)    
        mpp.legend(loc=0)

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)

    #end while

#end def Phisquare_plot


#------------------------------------------------------
# f(x,y)
# --------
# choice 21
#------------------------------------------------------
def deltaf2D_xy_plot(resu):

    """ delta f2D(x,y)"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.diag_time,give_default=2)

    idiag_x4 = np.int(resu.fdistribu_indx_diag[0][3])
    while (stime_diag != ''):

        delta_f_x1x2 = ( resu.fdistribu_x1x2[itime_diag,:,:] - 
                          resu.feq_xyvpar[idiag_x4,:,:] )

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = ax1.pcolormesh( resu.xgrid_2d[:,:], \
                                 resu.ygrid_2d[:,:], \
                                 delta_f_x1x2 ) ; 
        mpp.xlabel('x')
        mpp.ylabel('y')
        mpp.title('delta f(x,y)')
        fig.colorbar(p1)
        fig.show()

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)
    #end while

#end def deltaf2D_xy_plot


#------------------------------------------------------
# delta f(z,vpar)
# --------
# choice 22
#------------------------------------------------------
def deltaf2D_eta3vpar_plot(resu):

    """ f2D(eta3,vpar)"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.diag_time,give_default=2)

    Neta3  = np.size(resu.eta3_grid)
    Nvpar = np.size(resu.vpar_grid)
    feq_eta3vpar = np.zeros((Nvpar,Neta3))

    idiag_x1 = np.int(resu.fdistribu_indx_diag[0][0])
    idiag_x2 = np.int(resu.fdistribu_indx_diag[0][1])

    while (stime_diag != ''):

        for ivpar in np.r_[0:Nvpar]:
            for ieta3 in np.r_[0:Neta3]:
                feq_eta3vpar[ivpar,ieta3] = \
                    resu.feq_xyvpar[ivpar,idiag_x2,idiag_x1]
            #end for
        #end for

        delta_f_x3x4 = ( resu.fdistribu_x3x4[itime_diag,:,:] - 
                         feq_eta3vpar[:,:] )

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = ax1.pcolormesh(resu.eta3_grid[:],
                             resu.vpar_grid[:], 
                             delta_f_x3x4 ) 
        mpp.xlabel('eta3')
        mpp.ylabel('vpar')
        mpp.title('delta f(eta3,vpar)')
        fig.colorbar(p1)
        fig.show()

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)

    #end while

#end def deltaf2D_eta3vpar_plot


#------------------------------------------------------
# Norms
# --------
# choice 31
#------------------------------------------------------
def L1_L2_Linf_norms_plot(resu):

    """ L1-norm, L2-norm and Linf-norm"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.diag_time,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(6,6))
        #--> L1-norm
        ax1 = fig.add_subplot(1,3,1)
        p1  = mpp.plot(resu.diag_time,resu.L1_norm,'+-')
        mpp.xlabel('time')
        mpp.ylabel('L1-norm')
        mpp.grid(True)

        #--> L2-norm
        ax2 = fig.add_subplot(1,3,2)
        p2  = mpp.plot(resu.diag_time,resu.L2_norm,'+-')
        mpp.xlabel('time')
        mpp.ylabel('L2-norm')
        mpp.grid(True)

        #--> L3-norm
        ax3 = fig.add_subplot(1,3,3)
        p3  = mpp.plot(resu.diag_time,resu.Linf_norm,'+-')
        mpp.xlabel('time')
        mpp.ylabel('Linf-norm')
        mpp.grid(True)

        fig.show()

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)
    #end while
#end def L1_L2_Linf_norms_plot


#------------------------------------------------------
# Energy
# --------
# choice 32
#------------------------------------------------------
def energies_plot(resu):

    """ Kinetic and potential energies """
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.diag_time,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(6,6))
        ax1 = fig.add_subplot(1,1,1)
        mpp.hold(True)
        p1  = mpp.plot( resu.diag_time, 
                        resu.energy_kin,'-^g',
                        label='kinetic energy' )
        p2  = mpp.plot( resu.diag_time,
                        resu.energy_pot,'o-b',
                        label='potential energy' )
        #VG#p3  = mpp.plot( resu.diag_time,
        #VG#                resu.energy_tot,'-+m',
        #VG#                label='total energy')
        mpp.grid(True)
        mpp.legend()
        mpp.hold(False)

        fig.show()
        
        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)

    #end while

#end def energies_plot
