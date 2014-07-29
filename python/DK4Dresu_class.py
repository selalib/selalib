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

    #end def ___init__


    #---------------------------------------
    # Result reading
    #---------------------------------------
    def find_existing_results(self,ifirst_resufile=-1,ilast_resufile=-1):

        """ Result reading"""

        #--> Read all the existing results
        l_resufiles = glob.glob(self.Directory+'/DK4d_diag*.h5')
        l_resufiles = sorted(l_resufiles)
        self.fnames_resu = l_resufiles

        size_extent   = 3  # for '.h5'
        size_num      = 5  # size of the integer for diagnostic 
        l_diagnum     = map(lambda x : x[-(size_extent+size_num):-size_extent],l_resufiles)
        i_diagnum_min = int(min(l_diagnum)) 
        i_diagnum_max = int(max(l_diagnum)) 
        nb_diag       = i_diagnum_max-i_diagnum_min+1

        if ((ifirst_resufile==-1) or (ilast_resufile==-1)):
            [indx_first_ask,indx_last_ask] = \
                Ask_resultFiles(l_diagnum)
            self.ibeg  = indx_first_ask
            self.iend  = indx_last_ask
        else:
            self.ibeg = ifirst_resufile
            self.iend = ilast_resufile 
        #end if
        self.nfile = self.iend-self.ibeg+1

    #end def find_existing_results


    #--------------------------------------------------
    # Reading of the files 'conservation_laws.h5'
    #--------------------------------------------------
    def read_conservation_laws(self):

        """ Reading of the files 'conservation_laws.h5'"""

        #--> Load the file
        H = GYSut.loadHDF5(self.Prefix+'/conservation_laws.h5')
        for iname in H.keys:
            exec "self.%s = H.%s" % (iname,iname)
        #end for
    
    #end read_conservation_laws


    #--------------------------------------------------
    # Reading of the files 'DK4d_diag_d<num_diag>.h5'
    #--------------------------------------------------
    def read_diag(self):

        """ Reading of the files 'DK4d_diag_d<num_diag>.h5'"""

        #--> Find all the result files
        self.find_existing_results()

        #--> Load the first file
        H = GYSut.loadHDF5(self.fnames_resu[self.ibeg])
        for iname in H.keys:
            tupple_tmp = (self.nfile,)
            exec "sh = np.shape(H.%s)" % (iname)
            tupple_tmp = tupple_tmp + sh
            exec "self.%s = np.zeros(tupple_tmp)" % (iname)
    
        #--> Read all files
        ifile = 0
        for filename in self.fnames_resu[self.ibeg:self.iend+1]:
            filename_str = string.replace(filename,self.Prefix,'')
            print ' ---> Reading of the file : ' + self.Prefix + filename_str
            H = GYSut.loadHDF5(filename)
            for iname in H.keys:
                exec "self.%s[ifile,] = H.%s" % (iname,iname)
            #end for
            ifile = ifile + 1
        #end for
    #end def read_diag

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
    nb_DK4D = np.size(glob.glob(directory+'/DK4d_*.h5'))
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
    print ' *** GYSELA Diagnostics : choice of files *** '
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



