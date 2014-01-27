#-------------------------------------------------
#  file : DK4Ddiag_f2D.py
#  date : 2014-01-25
#   Diagnostics for DK4D code associated
#    to the saving of 2D electric potential
#-------------------------------------------------
"""
#  file : DK4Ddiag_f2D.py
#  date : 2014-01-27

Diagnostics for DK4D code associated
to the saving of distribution function
"""

import matplotlib as mp
import matplotlib.pyplot as mpp

import GYSutils
GYSut = GYSutils.impall('GYSutils')


#------------------------------------------------------
# f(x,y)
# --------
# choice 21
#------------------------------------------------------
def f2D_xy_plot(resu):

    """ f2D(x,y)"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.time_diag,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = ax1.pcolormesh(resu.xgrid_2d[itime_diag,:,:], \
                                 resu.ygrid_2d[itime_diag,:,:], \
                                 resu.f2d_xy[itime_diag,:,:]) 

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.time_diag)
    #end while
#end def f2D_xy_plot


#------------------------------------------------------
# f(z,vpar)
# --------
# choice 22
#------------------------------------------------------
def f2D_zvpar_plot(resu):

    """ f2D(z,vpar)"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu.time_diag,give_default=2)

    while (stime_diag != ''):

        fig = mpp.figure(figsize=(5,5))
        ax1 = fig.add_subplot(1,1,1)
        p1  = ax1.pcolormesh(resu.eta3_grid[itime_diag,:], \
                                 resu.vpar_grid[itime_diag,:], \
                                 resu.f2d_zvpar[itime_diag,:,:]) 

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu.time_diag)
    #end while

#end def f2D_zvpar_plot
