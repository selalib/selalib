#-------------------------------------------------
#  file : DK4Dbyb_main.py
#  date : 2014-08-25
#
# Main program for DK4D hybrid diagnostics
#-------------------------------------------------

""" Main program for DK4D hybrid diagnostics."""

import os
import GYSutils
GYSut = GYSutils.impall('DK4Dhyb_GYSutils')
DK4D  = GYSut.impall('DK4Dhyb_resu_class')
DK4Dd = GYSut.impall('DK4Dhyb_diag')

STOP = 0
Lresu_read = True

#--> Read DK4D results
if (Lresu_read):
    #--> Ask DK4D result directory
    DK4Ddir = DK4D.Ask_directory()
    print "DK4Ddir = " + str(DK4Ddir)
    if (DK4Ddir == -1):
        # There are NO DK4D results HDF5 files !
        STOP = 1
    #endif
#endif

if (STOP != 1):
    if (Lresu_read):
        resu = DK4D.DK4Dresu(DK4Ddir)
        resu.read_mesh()
        resu.read_equilibrium()
        resu.read_diags_files()
        resu.read_fdistribu_files()
        resu.read_Phi_files()
    #end if

    #--> DK4D diagnostics
    DK4Dd.DK4Ddiag(resu)
#end if

