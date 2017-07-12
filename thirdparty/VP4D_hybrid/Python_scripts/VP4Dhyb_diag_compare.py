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
GYSut = VP4Dhyb_GYSutils.impall('VP4Dhyb_GYSutils')


#-----------------------------------------------
# Diagnostic choices
#-----------------------------------------------
def VP4Ddiag(resu, nb_file, files_name):
    """ Diagnostic choices"""

    choice = '1';
    while choice != '' :
        print ' **** CONSERVATION LAWS ****'
        print '  *   1 : L1-Norm, L2-Norm and entropy'
        print '  *   2 : Energy'
        print ' '
        print ' **** LOG PLOT FOR CONVERGENCE ****'
        print '  *   11 : L1-Norm, L2-Norm and entropy'
        print '  *   12 : Energy'
        print ' ' 

        choice = raw_input(' choice ? : ');
        choice = GYSut.str2num(choice)

        if (choice != ''):
            VP4Ddiag_choice(choice, resu, nb_file, files_name)

    #end while

    print 'Bye bye!'

#end def VP4Ddiag


#-----------------------------------------------------------------
#  Execution of VP4D diagnostics according to the user choice
#-----------------------------------------------------------------
def VP4Ddiag_choice(ichoice, resu, nb_file, files_name):
    """
    Execution of GYSELA diagnostics according to the user choice
    """
    
    if (ichoice == 1):
        L1norm_L2norm_entropy_plot(resu, nb_file, files_name)
    elif (ichoice == 2):
        energies_plot(resu, nb_file, files_name)
    elif (ichoice == 11):
        L1norm_L2norm_entropy_log_plot(resu, nb_file, files_name)
    elif (ichoice == 12):
        energies_log_plot(resu, nb_file, files_name)

#end def VP4Ddiag_choice

#------------------------------------------------------
# Norms
# --------
# choice 1
#------------------------------------------------------
def L1norm_L2norm_entropy_plot(resu, nb_file, files_name):

    """ L1-norm, L2-norm and entropy"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu[0].diag_time,give_default=2)

    while (stime_diag != ''):
    
    	mpp.ion()

        mpp.figure(figsize=(6,6))
        mpp.hold(True)
        #--> L1-norm
        for i_file in range(nb_file):
            mpp.plot(resu[i_file].diag_time,resu[i_file].L1_norm,'+-', label=files_name[i_file, 1] + '; dt = ' + str(files_name[i_file, 2]) + '; nodes = ' + str(files_name[i_file,3]))
        #end for
        mpp.xlabel('time')
        mpp.ylabel('L1-norm')
        mpp.grid(True)

        mpp.hold(False)
        mpp.legend()

        mpp.figure(figsize=(6,6))
        mpp.hold(True)
        #--> L2-norm
        for i_file in range(nb_file):
            mpp.plot(resu[i_file].diag_time,resu[i_file].L2_norm,'+-', label=files_name[i_file, 1] + '; dt = ' + str(files_name[i_file, 2]) + '; nodes = ' + str(files_name[i_file,3]))
        #end for
        mpp.xlabel('time')
        mpp.ylabel('L2-norm')
        mpp.grid(True)

        mpp.hold(False)
        mpp.legend()

        mpp.figure(figsize=(6,6))
        mpp.hold(True)
        #--> entropy
        for i_file in range(nb_file):
            mpp.plot(resu[i_file].diag_time,resu[i_file].entropy_kin,'+-', label=files_name[i_file, 1] + '; dt = ' + str(files_name[i_file, 2]) + '; nodes = ' + str(files_name[i_file,3]))
        #end for
        mpp.xlabel('time')
        mpp.ylabel('kinetic entropy')
        mpp.grid(True)

        mpp.hold(False)
        mpp.legend()
	mpp.ioff()
#        mpp.show()

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu[0].diag_time)
    #end while
#end def L1norm_L2norm_entropy_plot


#------------------------------------------------------
# Energy
# --------
# choice 2
#------------------------------------------------------
def energies_plot(resu, nb_file, files_name):

    """ Total energy """
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu[0].diag_time,give_default=2)

    while (stime_diag != ''):
    
    	mpp.ion()

        fig = mpp.figure(figsize=(6,6))
        mpp.hold(True)

        for i_file in range(nb_file):

            Entot       = resu[i_file].energy_kin+resu[i_file].energy_pot
            Entot_0     = Entot-Entot[0]

            mpp.plot( resu[i_file].diag_time,
                        Entot_0,'-+',
                        label=files_name[i_file, 1] + '; dt = ' + str(files_name[i_file, 2]) + '; nodes = ' + str(files_name[i_file,3]))
        mpp.grid(True)
        mpp.xlabel('time')
        mpp.ylabel('Total energy')
        mpp.legend()
        mpp.hold(False)
	
	mpp.ioff()

#        fig.show()
        
        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu[0].diag_time)

    #end while

#end def energies_plot



#------------------------------------------------------
# Norms, log plot
# --------
# choice 11
#------------------------------------------------------
def L1norm_L2norm_entropy_log_plot(resu, nb_file, files_name):

    """ L1-norm, L2-norm and entropy"""
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu[0].diag_time,give_default=2)

    L1_norm = np.zeros(nb_file, dtype=float, order='C')
    L2_norm = np.zeros(nb_file, dtype=float, order='C')
    Entropy = np.zeros(nb_file, dtype=float, order='C')
    Time_step = np.zeros(nb_file, dtype=float, order='C')

    while (stime_diag != ''):

        for i_file in range(nb_file):

            L1_norm[i_file] = max(abs(resu[i_file].L1_norm[0:itime_diag+1]-resu[i_file].L1_norm[0]))[0]
            L2_norm[i_file] = max(abs(resu[i_file].L2_norm[0:itime_diag+1]-resu[i_file].L2_norm[0]))[0]
            Entropy[i_file] = max(abs(resu[i_file].entropy_kin[0:itime_diag+1]-resu[i_file].entropy_kin[0]))[0]
            Time_step[i_file] = files_name[i_file, 2]

        #end for

    	mpp.ion()

        mpp.figure(figsize=(6,6))
        mpp.hold(True)
        #--> L1-norm
        mpp.loglog(Time_step,L1_norm,'r+', label='Numerical values')
        
        poly_resu = np.polyfit( np.log10(Time_step), np.log10(L1_norm), 1 )
        
        str_legend = str(10**poly_resu[1]) + ' . dt^{' + str(poly_resu[0]) + '}'
        mpp.loglog(Time_step, 10**poly_resu[1]*Time_step**poly_resu[0], 'b-', label=str_legend)

        mpp.xlabel('time step dt')
        mpp.ylabel('L1-norm')
        mpp.grid(True)

        mpp.hold(False)
        mpp.legend()
        mpp.title('L1-norm maximum before t = ' + str(resu[0].diag_time[itime_diag][0]))

        mpp.figure(figsize=(6,6))
        mpp.hold(True)
        #--> L2-norm
        mpp.loglog(Time_step,L2_norm,'r+', label='Numerical values')
        
        poly_resu = np.polyfit( np.log10(Time_step), np.log10(L2_norm), 1 )

        str_legend_2 = str(10**poly_resu[1]) + ' . dt^{' + str(poly_resu[0]) + '}'
        mpp.loglog(Time_step, 10**poly_resu[1]*Time_step**poly_resu[0], 'b-', label=str_legend_2)

        mpp.xlabel('time')
        mpp.ylabel('L2-norm')
        mpp.grid(True)

        mpp.hold(False)
        mpp.legend()
        mpp.title('L2-norm maximum before t = ' + str(resu[0].diag_time[itime_diag][0]))

        mpp.figure(figsize=(6,6))
        mpp.hold(True)
        #--> entropy
        mpp.loglog(Time_step,Entropy,'r+', label='Numerical values')
        
        poly_resu = np.polyfit( np.log10(Time_step), np.log10(Entropy), 1 )

        str_legend_3 = str(10**poly_resu[1]) + ' . dt^{' + str(poly_resu[0]) + '}'
        mpp.loglog(Time_step, 10**poly_resu[1]*Time_step**poly_resu[0], 'b-', label=str_legend_3)

        mpp.xlabel('time')
        mpp.ylabel('kinetic entropy')
        mpp.grid(True)

        mpp.hold(False)
        mpp.legend()
        mpp.title('Entropy maximum before t = ' + str(resu[0].diag_time[itime_diag][0]))
	mpp.ioff()
#        mpp.show()

        print ' '
        print 'Time step : ', Time_step
        print 'L1-norm   : ', L1_norm
        print 'L2-norm   : ', L2_norm
        print 'Entropy   : ', Entropy
        print ' '

        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu[0].diag_time)
    #end while
#end def L1norm_L2norm_entropy_plot


#------------------------------------------------------
# Energy, log plot
# --------
# choice 12
#------------------------------------------------------
def energies_log_plot(resu, nb_file, files_name):

    """ Total energy """
    
    #--> Ask time
    [stime_diag,itime_diag] = \
        GYSut.Ask_time(resu[0].diag_time,give_default=2)

    while (stime_diag != ''):
    
    	mpp.ion()

        fig = mpp.figure(figsize=(6,6))
        mpp.hold(True)

        Entot_0 = np.zeros(nb_file, dtype=float,order='C')
        Time_step = np.zeros(nb_file, dtype=float,order='C')

        for i_file in range(nb_file):

            Entot             = resu[i_file].energy_kin+resu[i_file].energy_pot
            Entot_0[i_file]   = max(abs(Entot[0:itime_diag+1]-Entot[0]))[0]
            Time_step[i_file] = files_name[i_file, 2]

        #end for

        mpp.loglog( Time_step, Entot_0,'r+', label='Numerical values')

        poly_resu = np.polyfit( np.log10(Time_step), np.log10(Entot_0), 1 )

        str_legend_4 = str(10**poly_resu[1]) + ' . dt^{' + str(poly_resu[0]) + '}'
        mpp.loglog(Time_step, 10**poly_resu[1]*Time_step**poly_resu[0], 'b-', label=str_legend_4)

        mpp.grid(True)
        mpp.xlabel('time')
        mpp.ylabel('Total energy')
        mpp.legend()
        mpp.title('Total energy maximum before t = ' + str(resu[0].diag_time[itime_diag][0]))
        mpp.hold(False)
	
	mpp.ioff()

#        fig.show()

        print ' '
        print 'Time step    : ', Time_step
        print 'Total energy : ', Entot_0
        print ' '
        
        #--> Ask time        
        [stime_diag,itime_diag] = GYSut.Ask_time(resu[0].diag_time)

    #end while

#end def energies_plot
