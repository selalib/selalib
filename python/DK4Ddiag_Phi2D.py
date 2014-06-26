#-------------------------------------------------
#  file : DK4Ddiag_Phi2D.py
#  date : 2014-01-25
#   Diagnostics for DK4D code associated
#    to the saving of 2D electric potential
#-------------------------------------------------
"""
#  file : DK4Ddiag_Phi2D.py
#  date : 2014-01-25

Diagnostics for DK4D code associated
to the saving of 2D electric potential.
and the integral of the electric potentiel in square 
"""

import matplotlib as mp
import matplotlib.pyplot as mpp

import GYSutils
GYSut = GYSutils.impall('GYSutils')


#------------------------------------------------------
# Phi(x,y)
# --------
# choice 1
#------------------------------------------------------
def Phi2D_xy_plot(resu):
    """ Phi(x,y)"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.time_diag,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = ax1.pcolormesh(resu.xgrid_2d[itime_diag,:,:], \
                                 resu.ygrid_2d[itime_diag,:,:], \
                                 resu.phi2d_xy[itime_diag,:,:]) 

        fig.show()

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.time_diag)
    #end while
#end def Phi2D_xy_plot



#------------------------------------------------------
# log int Phi(eta1,eta2,eta3)**2 jac deta1 deta2 deta3
# --------
# choice 1
#------------------------------------------------------
def log_Int_phisquare_plot(resu):
    """ Phi(x,y)"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.time_diag,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = mpp.semilogy(resu.time_evol,resu.phisquare,'+-')
        mpp.xlabel('time')
        mpp.ylabel('intsquare')
        mpp.grid(True)
        fig.show()

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.time_diag)
    #end while
#end def Phi2D_xy_plot



