#-------------------------------------------------
#  file : VP4Dbyb_main.py
#  date : 2014-08-25
#
# Main program for VP4D hybrid diagnostics
#-------------------------------------------------

""" Main program for VP4D hybrid diagnostics."""

import os
import VP4Dhyb_GYSutils

#---- Define the Python script repository ----
s_cwd        = os.getcwd()
s_cwd_split  = s_cwd.split('VP4D_hybrid')
resudir_name = s_cwd.split('/')[-1]
python_path  = s_cwd_split[0] + 'VP4D_hybrid/Python_scripts'
home         = os.environ["HOME"]
python_path  = python_path.replace('/gpfsdata/vgrandgirard',home)
os_path      = os.sys.path
if (os_path.count(python_path)==0):
    os_path.insert(1,python_path)
#end if

#---- Load the useful modules ----
GYSut = VP4Dhyb_GYSutils.impall('VP4Dhyb_GYSutils')
VP4D  = GYSut.impall('VP4Dhyb_resu_class')
VP4Dd = GYSut.impall('VP4Dhyb_diag')
VP4Ds = GYSut.impall('VP4Dhyb_save')


#--> Read VP4D results
STOP = 0
Lresu_read = True
if (Lresu_read):
    #--> Ask VP4D result directory
    VP4Ddir = VP4D.Ask_directory()
    print "VP4Ddir = " + str(VP4Ddir)
    if (VP4Ddir == -1):
        # There are NO VP4D results HDF5 files !
        STOP = 1
    #endif
#endif

if (STOP != 1):
    if (Lresu_read):
        resu = VP4D.VP4Dresu(VP4Ddir)
        resu.read_mesh()
        resu.read_equilibrium()
        resu.read_diags_files()
        resu.read_fdistribu_files()
        resu.read_Phi_files()
    #end if

    #--> VP4D time evolution saving
    choice_save = raw_input('Do you want to save the time evolution data [default = n] ?: ')
    if ( (choice_save != '') & (choice_save =='y') ):
        FileNameResu = resudir_name+'_save'
        VP4Ds.SaveEvol_inHDF5(FileNameSave=FileNameResu,
                              CreateFile=True,
                              resu=resu,
                              save2D_evol=True)
    #end if

    #--> VP4D diagnostics
    VP4Dd.VP4Ddiag(resu)

#end if

