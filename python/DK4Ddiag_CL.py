#-------------------------------------------------
#  file : DK4Ddiag_CL.py
#  date : 2014-01-25
#   Diagnostics for DK4D code associated
#    to the conservation laws
#-------------------------------------------------
"""
#  file : DK4Ddiag_CL.py
#  date : 2014-01-27

Diagnostics for DK4D code associated
to the conservation laws
"""

import matplotlib as mp
import matplotlib.pyplot as mpp

import GYSutils
GYSut = GYSutils.impall('GYSutils')


#------------------------------------------------------
# Norms
# --------
# choice 31
#------------------------------------------------------
def L1_L2_Linf_norms_plot(resu):

    """ L1-norm, L2-norm and Linf-norm"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.time_diag,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(6,6))
        #--> L1-norm
        ax1 = fig.add_subplot(1,3,1)
        p1  = mpp.plot(resu.time_evol,resu.L1_norm,'+-')
        mpp.grid(True)

        #--> L2-norm
        ax2 = fig.add_subplot(1,3,2)
        p2  = mpp.plot(resu.time_evol,resu.L2_norm,'+-')
        mpp.xlabel('time')
        mpp.grid(True)

        #--> L3-norm
        ax3 = fig.add_subplot(1,3,3)
        p3  = mpp.plot(resu.time_evol,resu.Linf_norm,'+-')
        mpp.grid(True)

        fig.show()

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.time_diag)
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
        GYSut.Ask_time(resu.time_diag,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(6,6))
        ax1 = fig.add_subplot(1,1,1)
        mpp.hold(True)
        p1  = mpp.plot(resu.time_evol,resu.nrj_kin,'-^g',label='kinetic energy')
        p2  = mpp.plot(resu.time_evol,resu.nrj_pot,'o-b',label='potential energy')
        p3  = mpp.plot(resu.time_evol,resu.nrj_tot,'-+m',label='total energy')
        mpp.grid(True)
        mpp.legend()
        mpp.hold(False)

        fig.show()
        
        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.time_diag)

    #end while

#end def energies_plot
