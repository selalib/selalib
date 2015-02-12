#-------------------------------------------------
#  file : DK4Dmain.py
#  date : 01/02/2013
#
# Main program for DK4D diagnostics
#-------------------------------------------------

""" Main program for DK4D diagnostics."""

import os
import GYSutils
GYSut = GYSutils.impall('GYSutils')
DK4D  = GYSut.impall('DK4Dresu_class')
DK4Dd = GYSut.impall('DK4Ddiag_main')

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
        resu.read_conservation_laws()
        resu.read_diag()
    #end if

    #--> DK4D diagnostics
    DK4Dd.DK4Ddiag(resu)
#end if

