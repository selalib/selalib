"""
#  file : DK4Dresu_class.py
#  date : 2014-01-24

   Reading of all the result files
   and saving of the time evolution 
   of all DK4D outputs
"""

import glob
import os
import numpy as np
import string 

import GYSutils
GYSut   = GYSutils.impall('GYSutils')


#---------------------------------------
# GYSspecies class definition
#---------------------------------------
class DK4Dresu():    

    #---------------------------------------
    # DK4Dresu initialisation
    #---------------------------------------
    def __init__(self,DK4Ddir=''):

        #--> Initialization of the GYSELA result directory
        if (DK4Ddir==''):
            self.Directory = Ask_directory()
            if (status==-1):
                print '\033[0;31m ==> So, stop the script !!!!! \033[1;m '
                return
        else:
            self.Directory =  DK4Ddir
        print " "
        print " ********************************************************************"
        print " " + str(self.Directory)
        print " ********************************************************************"
        self.Prefix = self.Directory+'/'
        self.ibeg   = -1
        self.iend   = -1

    #end def ___init__


    #--------------------------------------------------
    # Reading of the file 'mesh.h5'
    #--------------------------------------------------
    def read_mesh(self):

        """ Reading of the file 'mesh.h5'"""

        #--> Load the file
        mesh_filename = self.Prefix+'/mesh.h5'
        print ' ---> Reading of the file : ' + mesh_filename
        H = GYSut.loadHDF5(mesh_filename)
        for iname in H.keys:
            exec "self.%s = H.%s" % (iname,iname)
        #end for

    #end def read_mesh


    #--------------------------------------------------
    # Reading of the file 'equilibrium.h5'
    #--------------------------------------------------
    def read_equilibrium(self):

        """ Reading of the file 'equilibrium.h5'"""

        #--> Load the file
        equilibrium_filename = self.Prefix+'/equilibrium.h5'
        print ' ---> Reading of the file : ' + equilibrium_filename 
        H = GYSut.loadHDF5(equilibrium_filename)
        for iname in H.keys:
            exec "self.%s = H.%s" % (iname,iname)
        #end for

    #end def read_equilibrium


    #--------------------------------------------------
    # Reading of the files 'diags_d<num_diag>.h5'
    #--------------------------------------------------
    def read_diags_files(self):

        """ Reading of the files 'diags_d<num_diag>.h5'"""

        #--> Find all the result files
        [fname_list_diags,nb_files,ibeg,iend] = \
            self.find_existing_results('diags_d*.h5')
        self.fname_list_diags = fname_list_diags
        self.nb_files = nb_files
        self.ibeg     = ibeg
        self.iend     = iend
        self.read_filename_list(fname_list_diags,nb_files,ibeg,iend)

    #end def read_diags_files


    #--------------------------------------------------
    # Reading of the files 'fdistribu_d<num_diag>.h5'
    #--------------------------------------------------
    def read_fdistribu_files(self):

        """ Reading of the files 'fdistribu_d<num_diag>.h5'"""

        #--> Find all the result files
        [fname_list_fdistribu,nb_files,ibeg,iend] = \
            self.find_existing_results('fdistribu_d*.h5', 
                                       ifirst_resufile=self.ibeg, 
                                       ilast_resufile=self.iend)

        self.fname_list_fdistribu = fname_list_fdistribu
        self.nb_files = nb_files
        self.ibeg     = ibeg
        self.iend     = iend
        self.read_filename_list(fname_list_fdistribu,nb_files,ibeg,iend)

    #end def read_fdistribu_files


    #--------------------------------------------------
    # Reading of the files 'Phi_d<num_diag>.h5'
    #--------------------------------------------------
    def read_Phi_files(self):

        """ Reading of the files 'Phi_d<num_diag>.h5'"""

        #--> Find all the result files
        [fname_list_Phi,nb_files,ibeg,iend] = \
            self.find_existing_results('Phi_d*.h5', 
                                       ifirst_resufile=self.ibeg, 
                                       ilast_resufile=self.iend)

        self.fname_list_Phi = fname_list_Phi
        self.nb_files = nb_files
        self.ibeg     = ibeg
        self.iend     = iend
        self.read_filename_list(fname_list_Phi,nb_files,ibeg,iend)

    #end def read_Phi_files


    #---------------------------------------
    # Result reading
    #---------------------------------------
    def find_existing_results(self, 
                              file_name_list,
                              ifirst_resufile=-1, 
                              ilast_resufile=-1):

        """ Result reading"""

        #--> Read all the existing results
        #VG#l_resufiles = glob.glob(self.Directory+'/diags*.h5')
        l_resufiles = glob.glob(self.Directory+'/'+file_name_list)
        l_resufiles = sorted(l_resufiles)

        size_extent   = 3  # for '.h5'
        size_num      = 5  # size of the integer for diagnostic 
        l_diagnum     = map(lambda x : x[-(size_extent+size_num):-size_extent],l_resufiles)
        i_diagnum_min = int(min(l_diagnum)) 
        i_diagnum_max = int(max(l_diagnum)) 
        nb_diag       = i_diagnum_max-i_diagnum_min+1

        if ((ifirst_resufile==-1) or (ilast_resufile==-1)):
            [indx_first_ask,indx_last_ask] = \
                Ask_resultFiles(l_diagnum)
            ibeg  = indx_first_ask
            iend  = indx_last_ask
        else:
            ibeg = ifirst_resufile
            iend = ilast_resufile 
        #end if
        nb_files = iend-ibeg+1

        return  [l_resufiles,nb_files,ibeg,iend]

    #end def find_existing_results


    #--------------------------------------------------
    # Create and fill all the arrays for the time 
    #  evolution of each quantities contained in the
    #  list of file names
    #--------------------------------------------------
    def read_filename_list(self,fname_list,nfile,ibeg,iend):

        """ """
        #--> Load the first file
        H = GYSut.loadHDF5(fname_list[ibeg])
        for iname in H.keys:
            tupple_tmp = (nfile,)
            exec "sh = np.shape(H.%s)" % (iname)
            tupple_tmp = tupple_tmp + sh
            exec "self.%s = np.zeros(tupple_tmp)" % (iname)
        
        #--> Read all files
        ifile = 0
        for filename in fname_list[ibeg:iend+1]:
            filename_str = string.replace(filename,self.Prefix,'')
            print ' ---> Reading of the file : ' + self.Prefix + filename_str
            H = GYSut.loadHDF5(filename)
            for iname in H.keys:
                exec "self.%s[ifile,] = H.%s" % (iname,iname)
            #end for
            ifile = ifile + 1
        #end for
        
    #end def read_filename_list

#end class DK4Dresu


#---------------------------------------------------------
# Ask which directory contains DK4D results to be 
#  analysed
#---------------------------------------------------------
def Ask_directory():

    """ Ask which directory contains DK4D results
    to be analysed"""

    directory = ''
    str_ask   = ' Directory for files reading'+ \
                '  [default = current] ? '
    directory = raw_input(str_ask);
    if (directory==''):
        directory = str(os.getcwd())
    nb_DK4D = np.size(glob.glob(directory+'/diags_*.h5'))
    if ( (nb_DK4D==0) ):
        print '\033[0;31m',\
            ' ==> There are NO DK4D results HDF5 files in :',\
            '\033[1;m',str(directory)
        status = -1
        return status
    return directory

#end def Ask_directory


#---------------------------------------------------------------------
# Ask which files will be treated in the GYSELA diagnostics
#  and create the list of file names required
#---------------------------------------------------------------------
def Ask_resultFiles(diagnum_list):

    """ Ask which files will be treated in the GYSELA diagnostics
    and create the list of file names required"""

    str_first_file = str(diagnum_list[0])
    str_last_file  = str(diagnum_list[-1])

    print '';
    print ' *** DK4D hybrid diagnostics : choice of files *** '
    str_first_ask = raw_input(' First file ? (between ' + \
                str_first_file + ' (default) and ' + \
                str_last_file + ' )')
    if (str_first_ask ==''): str_first_ask  = str_first_file
    str_last_ask  = raw_input(' Last file ?  (between ' + \
                 str_first_ask + ' and ' + str_last_file + \
                 ' (default) )')
    if (str_last_ask ==''): str_last_ask = str_last_file

    l_tmp          = filter(lambda x: int(x) >= int(str_first_ask), diagnum_list)
    indx_first_ask = diagnum_list.index(l_tmp[0])
    l_tmp          = filter(lambda x: int(x) <= int(str_last_ask), diagnum_list)
    indx_last_ask  = diagnum_list.index(l_tmp[-1])
    return [indx_first_ask,indx_last_ask]

#end def Ask_resultFiles



