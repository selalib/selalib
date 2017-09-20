#-------------------------------------------------
#  file : VP4Dhyb_diag.py
#  date : 2014-08-25
#
#   Diagnostics for VP4D code
#-------------------------------------------------
import VP4Dhyb_GYSutils
import matplotlib.pyplot as mpp
import numpy as np
import scipy as sp
import GR_Freq_Computation as GRF
GYSut = VP4Dhyb_GYSutils.impall('VP4Dhyb_GYSutils')

#-----------------------------------------------
# Diagnostic choices
#-----------------------------------------------
def VP4Ddiag(resu):
    """ Diagnostic choices"""

    choice = '1';
    while choice != '' :
        print ' **** MESH **** '
        print '  *   0 : mesh '
        print ' '
        print ' **** EQUILIBRIUM ****'
        print '  *   1 : equilibrium distribution function'
        print ' '
        print ' **** PHI2D DIAGNOSTICS ****'
        print '  *   11  : Phi(x,y)'
        print '  *   12  : int Phi(x,y,z)**2 jac dx dy dz'
        print '  *   13  : int E**2 jac dx dy dz'
        print '  *   112 : int Phi(x,y,z)**2 jac dx dy dz (new version)' 
        print '  *   113 : int E**2 jac dx dy dz (new version)'
        print ' '

        print ' **** F2D DIAGNOSTICS ****'
        print '  *   21 : f2D(x,y)'
        print '  *   22 : delta f2D(x,y)'
        print '  *   23 : f2D(vx,vy)'
        print '  *   24 : delta f2D(vx,vy)'
        print ' '

        print ' **** CONSERVATION LAWS ****'
        print '  *   31 : L1-Norm, L2-Norm and entropy'
        print '  *   32 : Energy'
        print ' '

        choice = raw_input(' choice ? : ');
        choice = GYSut.str2num(choice)

        if (choice != ''):
            VP4Ddiag_choice(choice,resu)

    #end while

    print 'Bye bye!'

#end def VP4Ddiag


#-----------------------------------------------------------------
#  Execution of VP4D diagnostics according to the user choice
#-----------------------------------------------------------------
def VP4Ddiag_choice(ichoice,resu):
    """
    Execution of GYSELA diagnostics according to the user choice
    """
    
    if (ichoice == 0):
        mesh_plot(resu)
    elif (ichoice == 1):
        equilibrium_fdistribu_plot(resu)
    elif (ichoice == 11):
        Phi2D_xy_plot(resu)
    elif (ichoice == 12):
        Phisquare_plot(resu)
    elif (ichoice == 13):
        Esquare_plot(resu)
    elif (ichoice == 112):
        Phisquare_plot_GFerriere(resu)
    elif (ichoice == 113):
        Esquare_plot_GFerriere(resu)
    elif (ichoice == 21):
        f2D_xy_plot(resu)
    elif (ichoice == 22):
        deltaf2D_xy_plot(resu)
    elif (ichoice == 23):
        f2D_vxvy_plot(resu)
    elif (ichoice == 24):
        deltaf2D_vxvy_plot(resu)
    elif (ichoice == 31):
        L1norm_L2norm_entropy_plot(resu)
    elif (ichoice == 32):
        energies_plot(resu)

#end def VP4Ddiag_choice


#------------------------------------------------------
# Plot the mesh
# --------
# choice 0
#------------------------------------------------------
def mesh_plot(resu):
    """ Plot the physical mesh associated to the simulation """

    fig = mpp.figure(figsize=(6,6))
    mpp.hold(True)
    Nx = np.shape(resu.xgrid_2d)[0]
    Ny = np.shape(resu.xgrid_2d)[1]
    for i in np.r_[0:Nx]:
        mpp.plot(resu.xgrid_2d[i,:],resu.ygrid_2d[i,:],'k')
    #end for
    for j in np.r_[0:Ny]:
        mpp.plot(resu.xgrid_2d[:,j],resu.ygrid_2d[:,j],'k')
    #end for
    mpp.hold(False)
    fig.show()

#end def mesh_plot
#------------------------------------------------------


#------------------------------------------------------
# Plot the equilibrium distribution function
# --------
# choice 1
#------------------------------------------------------
def equilibrium_fdistribu_plot(resu):

    """ Equilibrium distribution function """

    fig = mpp.figure(figsize=(6,6))

    ax1 = fig.add_subplot(1,1,1)
    p1 = ax1.pcolor( 
        resu.vx_grid,
        resu.vy_grid,
        resu.feq_vxvy[:,:] )
    fig.colorbar(p1)
    mpp.xlabel('vx')
    mpp.ylabel('vy')
    mpp.title('feq(vx,vy)')
    mpp.grid(True)
    fig.show()

    fig.show()

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
    """ \int Phi(x,y)**2 dx dy"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.diag_time,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = mpp.semilogy( 
            resu.diag_time[0:itime_diag + 1], 
            np.sqrt(resu.Phisquare[0:itime_diag + 1]),'+-' )
        mpp.xlabel('time')
        mpp.ylabel('int Phi^2')
        mpp.grid(True)
        fig.show()

        #--> Computation of the frequence and the growth rate
        print ' ##### parameters for the linear interpolation and for the precision of the Fourier transform for omega #####'
        [stime_init_interp,itime_init_interp] = GYSut.Ask_time(resu.diag_time, \
               ask_msg='Time for interpolation beginning',give_default=1)
        [stime_end_interp,itime_end_interp] = GYSut.Ask_time(resu.diag_time, \
               first_indx_time=itime_init_interp,               
               ask_msg='Time for interpolation end',give_default=2)
        d_itime_interp  = float(raw_input('Precision of the time for interpolation beginning and end : '))
        Nk              = int(raw_input('Number of points in frequence for the Fourier transform : '))
        dt              = resu.diag_time[1] - resu.diag_time[0]
        imax            = int(d_itime_interp / dt)
        time_tmp        = resu.diag_time[itime_init_interp:itime_end_interp+1].T[0]
        Phi00           = resu.Phi_x1x2[itime_init_interp:itime_end_interp+1,0,0]
        [TFFPhi, Ffreq] = GYSut.Fourier1D( Phi00, time_tmp )
        omegaapp        = np.abs(Ffreq[np.max(np.nonzero(np.abs(TFFPhi)==np.max(np.abs(TFFPhi))))])
        Results         = np.zeros((2*imax+1, 2*imax+1, 5), dtype=float, order='C')
        for i in range(2 * imax + 1):
            for j in range(2 * imax + 1):
                itime_init_interp0 = itime_init_interp + (i - imax)
                itime_end_interp0 = itime_end_interp + (j - imax)
                time_tmp        = resu.diag_time[itime_init_interp0:itime_end_interp0+1].T[0]
                Phi00           = resu.Phi_x1x2[itime_init_interp0:itime_end_interp0+1,0,0]
                Lx              = time_tmp[-1] - time_tmp[0]
                [TFPhi, freq]   = GYSut.TFourier1D( Phi00, time_tmp, omegaapp + 2 * sp.pi / Lx, max([1, omegaapp - 2 * sp.pi / Lx]), Nk)
                omega           = np.abs(freq[np.max(np.nonzero(np.abs(TFPhi)==np.max(np.abs(TFPhi))))])
                nend            = int(sp.pi / (omega * dt))
                Phisquare_ct    = np.sqrt(resu.Phisquare[itime_init_interp0 + range(0,nend+1)].T[0])
                iphase          = np.nonzero(Phisquare_ct==np.min(Phisquare_ct))[0]
                # Computation of the phase
                #--> approximied phase
                phase           = time_tmp[iphase] + dt/100
                #--> modified phase
                #y1              = Phi00[iphase]
                #y2              = Phi00[iphase+1]
                #y3              = Phi00[iphase-1]
                #phase           = time_tmp[iphase] + 2*(y2+y3-2*y1)/(omega**2*dt*(y2-y3))
                # end of computation of the phase
                #print phase
                cos_omega_t     = np.abs(np.sin(omega*(time_tmp-phase)))
                Phisquare_tmp   = np.sqrt(resu.Phisquare[itime_init_interp0:itime_end_interp0+1].T[0])/cos_omega_t
                poly_resu       = np.polyfit( time_tmp, np.log(Phisquare_tmp), 1, full='True' )
                [p,s]           = poly_resu[0]
                residuals       = poly_resu[1]
                rank            = poly_resu[2] 
                singular_values = poly_resu[3]
                rcond           = poly_resu[4]
                Results[i,j,0]  = np.linalg.norm(np.exp(s)*np.exp(p*time_tmp)*np.abs(np.sin(omega*(time_tmp-phase))) - Phisquare_tmp)/np.sqrt(Lx)
                Results[i,j,1:5] = [omega, p, s, residuals]
            #end for
            
        #end for
        L2_diff         = np.amin(Results[:,:,0])
        Ind_best        = np.nonzero(Results[:,:,0] == L2_diff)
        ib              = Ind_best[0][0]
        jb              = Ind_best[1][0]
        omega           = Results[ib,jb,1]
        p               = Results[ib,jb,2]
        s               = Results[ib,jb,3]
        residuals       = Results[ib,jb,4]
        print '=> frequence              = ', omega
        print '=> growth_rate            = ', p
        print '=> residuals              = ', residuals
        print '=> Difference in L^2-norm = ', L2_diff

        mpp.hold(True)
        str_legend = r'$\gamma = $' + str(p) + '\n residual = ' + str(residuals)
        str_legend2 = str(np.exp(s)) + '. exp(' + str(p) + ' t) . cos(' + str(omega) + ' t + phase)'
        mpp.semilogy( resu.diag_time[0:itime_end_interp+1], 
                      np.exp(s)*np.exp(p*resu.diag_time[0:itime_end_interp+1]),'r',
                      label=str_legend )

        mpp.semilogy( resu.diag_time[0:itime_end_interp+1],
                      np.exp(s)*np.exp(p*resu.diag_time[0:itime_end_interp+1])*np.abs(np.sin(omega*(resu.diag_time[0:itime_end_interp+1]-phase)+dt/100)),'g',
                      label=str_legend2 )
        ax1.set_ylim(np.exp(s+p* resu.diag_time[-1].T[0])/10, 10)
        mpp.hold(False)    
        mpp.legend(loc=0)

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)

    #end while

#end def Phisquare_plot

#------------------------------------------------------
# log int Phi(eta1,eta2,eta3)**2 jac deta1 deta2 deta3
# --------
# choice 12 (developed by Guillaume Ferriere 06/2015)
#------------------------------------------------------
def Phisquare_plot_GFerriere(resu):
    """ \int Phi(x,y)**2 dx dy"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.diag_time,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = mpp.semilogy( 
            resu.diag_time[0:itime_diag + 1], 
            np.sqrt(resu.Phisquare[0:itime_diag + 1]),'+-' )
        mpp.xlabel('time')
        mpp.ylabel('int Phi^2')
        mpp.grid(True)
        fig.show()

        #--> Computation of the frequence and the growth rate
        print ' ##### parameters for the linear interpolation and for the precision of the Fourier transform for omega #####'
        [stime_init_interp,itime_init_interp] = GYSut.Ask_time(resu.diag_time, \
               ask_msg='Time for interpolation beginning',give_default=1)
        [stime_end_interp,itime_end_interp] = GYSut.Ask_time(resu.diag_time, \
               first_indx_time=itime_init_interp,               
               ask_msg='Time for interpolation end',give_default=2)
        #'Number of points in frequence for the Fourier transform : '))
        Nk              = int(raw_input('Number of points for omega : '))
        Phi00           = resu.Phi_x1x2[:,0,0]
        Phisquare       = resu.Phisquare
        diag_time       = resu.diag_time
        [omega, p, s, L2_diff, phase, residuals] = GRF.GR_Freq_Comp(Phi00, Phisquare, itime_init_interp, itime_end_interp, Nk, diag_time, 0, 0)
        print '=> frequence              = ', omega
        print '=> growth_rate            = ', p
        print '=> residuals              = ', residuals
        print '=> Difference in L^2-norm = ', L2_diff

        mpp.hold(True)
        str_legend = r'$\gamma = $' + str(p) + '\n residual = ' + str(residuals)
        str_legend2 = str(np.exp(s)) + '. exp(' + str(p) + ' t) . cos(' + str(omega) + ' t + phase)'
        mpp.semilogy( resu.diag_time[0:itime_end_interp+1], 
                      np.exp(s)*np.exp(p*resu.diag_time[0:itime_end_interp+1]),'r',
                      label=str_legend )

        mpp.semilogy( resu.diag_time[0:itime_end_interp+1],
                      np.exp(s)*np.exp(p*resu.diag_time[0:itime_end_interp+1])*np.abs(np.sin(omega*(resu.diag_time[0:itime_end_interp+1]-phase))),'g',
                      label=str_legend2 )
        ax1.set_ylim(np.exp(s+p* resu.diag_time[-1].T[0])/10, 10)
        mpp.hold(False)    
        mpp.legend(loc=0)

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)

    #end while

#end def Phisquare_plot_GFerriere


#------------------------------------------------------
# log int E(eta1,eta2,eta3)**2 jac deta1 deta2 deta3
# --------
# choice 13
#------------------------------------------------------
def Esquare_plot(resu):
    """ \int E(x,y)^2 dx dy"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.diag_time,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = mpp.semilogy( 
            resu.diag_time[0:itime_diag + 1], 
            np.sqrt(resu.energy_pot[0:itime_diag + 1]),'+-' )
        mpp.xlabel('time')
        mpp.ylabel('sqrt(int E^2 dx dy)')
        mpp.grid(True)
        fig.show()

        #--> Computation of the growth rate
        print ' ##### parameters for the linear interpolation and for the precision of the Fourier transform for omega #####'
        [stime_init_interp,itime_init_interp] = GYSut.Ask_time(resu.diag_time, \
               ask_msg='Time for interpolation beginning',give_default=1)
        [stime_end_interp,itime_end_interp] = GYSut.Ask_time(resu.diag_time, \
               first_indx_time=itime_init_interp,               
               ask_msg='Time for interpolation end',give_default=2)
        Nk = int(raw_input('Number of points in frequence for the Fourier transform : '))
        dt              = time_tmp[1] - time_tmp[0]
        time_tmp        = resu.diag_time[itime_init_interp:itime_end_interp+1].T[0]
        Energyp         = resu.energy_pot[itime_init_interp:itime_end_interp+1].T[0]
        [TFFE, Efreq]   = GYSut.Fourier1D( Energyp, time_tmp )
        omegaapp        = np.abs(Efreq[np.max(np.nonzero(np.abs(TFFE)*(np.abs(Efreq)>1)==np.max(np.abs(TFFE)*(np.abs(Efreq)>1))))])
        Lx              = time_tmp[-1] - time_tmp[0]
        [TFE, freq]     = GYSut.TFourier1D( Energyp, time_tmp, omegaapp + 2 * sp.pi / Lx, max([1, omegaapp - 2 * sp.pi / Lx]), Nk)
        omega           = np.abs(freq[np.max(np.nonzero(np.abs(TFE)==np.max(np.abs(TFE))))])/2
        nend            = int(sp.pi / (omega * dt))
        Esquare_ct      = np.sqrt(Energyp[0:nend+1])
        phase           = time_tmp[np.nonzero(Esquare_ct==np.min(Esquare_ct))]
        sqrt_Epot_tmp   = np.sqrt(Energyp)/np.abs(np.sin(omega*(time_tmp-phase)+dt/100))
        poly_resu       = np.polyfit( time_tmp, np.log(sqrt_Epot_tmp), 1, full='True' )
        [p,s]           = poly_resu[0]
        residuals       = poly_resu[1]
        rank            = poly_resu[2] 
        singular_values = poly_resu[3]
        rcond           = poly_resu[4]
        print '=> frequence = ', omega
        print '=> growth_rate = ', p
        print '=> residuals   = ', residuals

        mpp.hold(True)
        str_legend = r'$\gamma = $' + str(p) + '\n residual = ' + str(residuals)
        str_legend2 = str(np.exp(s)) + ' . exp(' + str(p) + ' t) . cos(' + str(omega) + ' t + phase)'
        mpp.semilogy( resu.diag_time[0:itime_end_interp+1], 
                      np.exp(s)*np.exp(p*resu.diag_time[0:itime_end_interp+1]),'r',
                      label=str_legend )
        
        mpp.semilogy( resu.diag_time[0:itime_end_interp+1],
                      np.exp(s)*np.exp(p*resu.diag_time[0:itime_end_interp+1].T[0])*np.abs(np.sin(omega*(resu.diag_time[0:itime_end_interp+1].T[0] -phase)+dt/100)),'g',
                      label=str_legend2 )
        
        ax1.set_ylim(np.exp(s + p* resu.diag_time[-1].T[0])/10, 10)
        mpp.hold(False)    
        mpp.legend(loc=0)
        
        fig2 = mpp.figure(figsize=(5,5))
        ax2 = fig.add_subplot(1,1,1)
        p2 = mpp.plot(freq, np.abs(TFE), 'b')
        mpp.title('FFT of \int E**2 dx dy')
        mpp.grid(True)
        fig2.show()


        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)

    #end while

#end def Esquare_plot


#------------------------------------------------------
# log int E(eta1,eta2,eta3)**2 jac deta1 deta2 deta3
# --------
# choice 13 (developed by Guillaume Ferriere 06/2015)
#------------------------------------------------------
def Esquare_plot_GFerriere(resu):
    """ \int E(x,y)^2 dx dy"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.diag_time,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = mpp.semilogy( 
            resu.diag_time[0:itime_diag + 1], 
            np.sqrt(resu.energy_pot[0:itime_diag + 1]),'+-' )
        mpp.xlabel('time')
        mpp.ylabel('sqrt(int E^2 dx dy)')
        mpp.grid(True)
        fig.show()

        #--> Computation of the growth rate
        print ' ##### parameters for the linear interpolation and for the precision of the Fourier transform for omega #####'
        [stime_init_interp,itime_init_interp] = GYSut.Ask_time(resu.diag_time, \
               ask_msg='Time for interpolation beginning',give_default=1)
        [stime_end_interp,itime_end_interp] = GYSut.Ask_time(resu.diag_time, \
               first_indx_time=itime_init_interp,               
               ask_msg='Time for interpolation end',give_default=2)
#        d_itime_interp = float(raw_input('Precision of the time for interpolation beginning and end : '))
        Nk = int(raw_input('Number of points for omega : '))                                             #'Number of points in frequence for the Fourier transform : '))
        energypot   = resu.energy_pot
        [omega, p, s, L2_diff, phase, residuals] = GRF.GR_Freq_Comp(energypot.T[0], energypot, itime_init_interp, itime_end_interp, Nk, resu.diag_time, 1, 1)
        print '=> frequence              = ', omega
        print '=> growth_rate            = ', p
        print '=> residuals              = ', residuals
        print '=> Difference in L^2-norm = ', L2_diff
        mpp.hold(True)
        str_legend = r'$\gamma = $' + str(p) + '\n residual = ' + str(residuals)
        str_legend2 = str(np.exp(s)) + ' . exp(' + str(p) + ' t) . cos(' + str(omega) + ' t + phase)'
        mpp.semilogy( resu.diag_time[0:itime_end_interp+1], 
                      np.exp(s)*np.exp(p*resu.diag_time[0:itime_end_interp+1]),'r',
                      label=str_legend )
        
        mpp.semilogy( resu.diag_time[0:itime_end_interp+1],
                      np.exp(s)*np.exp(p*resu.diag_time[0:itime_end_interp+1].T[0])*np.abs(np.sin(omega*(resu.diag_time[0:itime_end_interp+1].T[0] -phase))),'g',
                      label=str_legend2 )
        
        ax1.set_ylim(np.exp(s + p* resu.diag_time[-1].T[0])/10, 10)
        mpp.hold(False)    
        mpp.legend(loc=0)


        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)

    #end while

#end def Esquare_plot_GFerriere


#------------------------------------------------------
# f(x,y)
# --------
# choice 21
#------------------------------------------------------
def f2D_xy_plot(resu):

    """ f2D(x,y)"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.diag_time,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = ax1.pcolormesh( resu.xgrid_2d[:,:], 
                              resu.ygrid_2d[:,:], 
                              resu.fdistribu_x1x2[itime_diag,:,:] ) ; 
        mpp.xlabel('x')
        mpp.ylabel('y')
        mpp.title('f(x,y)')
        fig.colorbar(p1)
        fig.show()

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)
    #end while

#end def f2D_xy_plot


#------------------------------------------------------
# delta f(x,y)
# --------
# choice 22
#------------------------------------------------------
def deltaf2D_xy_plot(resu):

    """ delta f2D(x,y)"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.diag_time,give_default=2)

    idiag_x3 = np.int(resu.fdistribu_indx_diag[0][2])
    idiag_x4 = np.int(resu.fdistribu_indx_diag[0][3])
    while (stime_diag != ''):

        delta_f_x1x2 = ( resu.fdistribu_x1x2[itime_diag,:,:] - 
                          resu.feq_vxvy[idiag_x3,idiag_x4] )

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
# f(vx,vy)
# --------
# choice 23
#------------------------------------------------------
def f2D_vxvy_plot(resu):

    """ f2D(vx,vy)"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.diag_time,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = ax1.pcolormesh(resu.vx_grid[:],
                             resu.vy_grid[:], 
                             resu.fdistribu_x3x4[itime_diag,:,:] ) 
        mpp.xlabel('vx')
        mpp.ylabel('vy')
        mpp.title('f(vx,vy)')
        fig.colorbar(p1)
        fig.show()

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)

    #end while

#end def f2D_vxvy_plot


#------------------------------------------------------
# delta f(vx,vy)
# --------
# choice 24
#------------------------------------------------------
def deltaf2D_vxvy_plot(resu):

    """ delta f2D(vx,vy)"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.diag_time,give_default=2)

    idiag_x1 = np.int(resu.fdistribu_indx_diag[0][0])
    idiag_x2 = np.int(resu.fdistribu_indx_diag[0][1])

    while (stime_diag != ''):

        delta_f_x3x4 = ( resu.fdistribu_x3x4[itime_diag,:,:] - 
                         resu.feq_vxvy[:,:] )

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = ax1.pcolormesh(resu.vx_grid[:],
                             resu.vy_grid[:], 
                             delta_f_x3x4 ) 
        mpp.xlabel('vx')
        mpp.ylabel('vy')
        mpp.title('delta f(vx,vy)')
        fig.colorbar(p1)
        fig.show()

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)

    #end while

#end def deltaf2D_vxvy_plot


#------------------------------------------------------
# Norms
# --------
# choice 31
#------------------------------------------------------
def L1norm_L2norm_entropy_plot(resu):

    """ L1-norm, L2-norm and entropy"""
    
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

        #--> entropy
        ax3 = fig.add_subplot(1,3,3)
        p3  = mpp.plot(resu.diag_time,resu.entropy_kin,'+-')
        mpp.xlabel('time')
        mpp.ylabel('kinetic entropy')
        mpp.grid(True)

        fig.show()

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)
    #end while
#end def L1norm_L2norm_entropy_plot


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

        Enkin = resu.energy_kin - resu.energy_kin[0]
        delta_Enkin = Enkin/resu.energy_kin[0]
        Enpot = resu.energy_pot - resu.energy_pot[0]
        delta_Enpot = Enpot/resu.energy_pot[0]
        Entot       = resu.energy_kin+resu.energy_pot
        Entot_0     = Entot-Entot[0]
        delta_Entot = Entot_0/Entot[0]

        fig = mpp.figure(figsize=(6,6))
        ax1 = fig.add_subplot(1,1,1)
        mpp.hold(True)
        p1  = mpp.plot( resu.diag_time, 
                        delta_Enkin,'-^g',
                        label='kinetic energy' )
        p2  = mpp.plot( resu.diag_time,
                        delta_Enpot,'o-b',
                        label='potential energy' )
        p3  = mpp.plot( resu.diag_time,
                        delta_Entot,'-+m',
                        label='total energy')
        mpp.grid(True)
        mpp.legend()
        mpp.title("Relative difference")
        mpp.hold(False)

        fig.show()

        fig2 = mpp.figure(figsize=(6,6))
        ax2 = fig2.add_subplot(3,1,1)
        mpp.hold(True)
        p12  = mpp.plot( resu.diag_time, 
                        resu.energy_kin,'-^g',
                        label='kinetic energy' )
        mpp.grid(True)
        mpp.legend()
        ax3 = fig2.add_subplot(3,1,2)
        p22  = mpp.plot( resu.diag_time,
                        resu.energy_pot,'o-b',
                        label='potential energy' )
        mpp.grid(True)
        mpp.legend()
        ax4 = fig2.add_subplot(3,1,3)
        p32  = mpp.plot( resu.diag_time,
                        Entot,'-+m',
                        label='total energy')
        mpp.grid(True)
        mpp.legend()
        mpp.hold(False)

        fig2.show()
        
        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.diag_time)

    #end while

#end def energies_plot
