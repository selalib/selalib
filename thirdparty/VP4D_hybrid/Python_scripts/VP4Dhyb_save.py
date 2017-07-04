#-------------------------------------------------
#  file : VP4Dbyb_save.py
#  date : 2015-09-06
#
# Routine used to save results in HDF5 format
#-------------------------------------------------
import glob
import h5py as h5
import os
import sys

#---------------------------------------------
#
#---------------------------------------------
def SaveEvol_inHDF5(FileNameSave,CreateFile,resu,save2D_evol):

    HDF5filename_2d = FileNameSave + "_2d.h5"
    if ( (CreateFile) & os.path.exists(HDF5filename_2d)):
        os.system('rm '+HDF5filename_2d)
    #end if
    h5f_2d = h5.File(HDF5filename_2d,'a')

    if (save2D_evol):
        HDF5filename_3d = FileNameSave + "_3d.h5"
        if ( (CreateFile) & os.path.exists(HDF5filename_3d)):
            os.system('rm '+HDF5filename_3d)
        #end if
        h5f_3d = h5.File(HDF5filename_3d,'a')
    #end if

    dic_resu = dir(resu)
    for varname in dic_resu:
        exec "var = resu.%s" %(varname)
        var_type = type(var)
        if ( 'ndim' in dir(var) ):
            var_dim = var.ndim
            if ( var_dim <= 2 ):
                exec "h5f_2d.create_dataset('%s',data=var) " %(varname)
            else:
                if (save2D_evol):
                    exec "h5f_3d.create_dataset('%s',data=var) " %(varname)
                #end if
        #end if
    #end for

    h5f_2d.close()
    print " "
    print "==> Saving of evolution of 0 to 1D quantities in " + HDF5filename_2d
    if (save2D_evol):
        h5f_3d.close()
        print "==> Saving of evolution of 2D quantities in " + HDF5filename_3d
    print " "

#end def SaveEvol_inHDF5

