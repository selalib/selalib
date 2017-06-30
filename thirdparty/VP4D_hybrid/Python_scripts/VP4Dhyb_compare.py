#-------------------------------------------------
#  file : VP4Dbyb_main.py
#  date : 2014-08-25
#
# Main program for VP4D hybrid diagnostics
#-------------------------------------------------

""" Main program for VP4D hybrid diagnostics."""

import os
import VP4Dhyb_GYSutils
GYSut = VP4Dhyb_GYSutils.impall('VP4Dhyb_GYSutils')
VP4D  = GYSut.impall('VP4Dhyb_resu_class')
VP4Ddc = GYSut.impall('VP4Dhyb_diag_compare')

STOP = 0
Lresu_read = True

resu = {}
files_name = {}
i_file = 0

while (STOP != 1):
    #--> Read VP4D results
    if (Lresu_read):
        #--> Ask VP4D result directory
        VP4Ddir = VP4D.Ask_directory()
        print "VP4Ddir = " + str(VP4Ddir)
        if (VP4Ddir == -1):
            # There are NO VP4D results HDF5 files !
            STOP = 1
        #endif
    #endif

    files_name[i_file, 1] = raw_input('Which case is it ? (identity/Collela/...) : ')
    files_name[i_file, 2] = float(raw_input('Time step ? : '))
    files_name[i_file, 3] = int(raw_input('Number of nodes ? : '))

    if (STOP != 1):
        if (Lresu_read):
            resu[i_file] = VP4D.VP4Dresu(VP4Ddir)
            resu[i_file].read_mesh()
            resu[i_file].read_equilibrium()
            resu[i_file].read_diags_files()
            resu[i_file].read_fdistribu_files()
            resu[i_file].read_Phi_files()
        #end if

        answer = raw_input('Is there another file to be analysed ? (answer y or n) : ')
    
        if answer == 'y':
            i_file = i_file + 1
        else:
            STOP = 1
        #endif
    #endif
#end if

nb_file = i_file + 1


#--> VP4D diagnostics comparison
if (VP4Ddir != 1):
    VP4Ddc.VP4Ddiag(resu, nb_file, files_name)
#endif
