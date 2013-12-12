#! /usr/bin/env python

#********************************************************
# File: setup.py
# Date: 2013-12-11
#
# Define the compiling procedure for PYRELDISP package
#   which can be installed with the command:
#     python setup.py install
#
# Rk: The compiler options are defined in the file
#       setup.cfg
#********************************************************
from numpy.distutils.core import setup, Extension
import sys
import copy
import os

extra_link_args=[]
if sys.platform=='darwin':
    extra_link_args=['-framework','veclib']


#----------------------------------------------------------------------
#---> Define the list of files which must be compiled in 
#      the repository zeal_dir containing the ZEAL package
#----------------------------------------------------------------------
zeal_dir = 'zeal_src/'
zeal_files_list1=['prec.f90','error.f90', 'zeal_input.f90']
zeal_files_list2=['math_tools.f90','integr_input.f90', 'quad.f90',
                  'zeros.f90','split.f90','refine.f90', 
                  'zeal4py.f90','dqag.f', 'dqagx.f', 'wofz.f']
zeal_files_list1 = map(lambda x: zeal_dir+x,zeal_files_list1)
zeal_files_list2 = map(lambda x: zeal_dir+x,zeal_files_list2)


#----------------------------------------------------------------------
#---> Define the list of files dedicated to the
#      calculation of zeros of dispersion relation
#      and fill a dictionary containing the different
#      cases which could be treated (e.g landau, bump-on-tail, ...)
#
# Rk: To be able to create a new library you need to add the
#       name of your new case (let's call it 'newcase')
#       in the following directionary 'pyreldisp_dic' with 
#       a list of the associated files which must be compiled.
#      Ex: pyreldisp_dic['newcase'] = ['fun_newcase.f90']  
#----------------------------------------------------------------------
reldisp_files_list = 'PlasmaFunctions.f90'
pyreldisp_dic={}
pyreldisp_dic['landau']  = [reldisp_files_list,'fun_landau.f90']
pyreldisp_dic['BOT']     = [reldisp_files_list,'fun_BOT.f90']
pyreldisp_dic['gyrokin'] = [reldisp_files_list,'fun_gyrokin.f90']


#----------------------------------------------------------------------
#---> Create the names of the different python packages which
#      will be constructed and define the list of 
#      source files associated
# Rk: One python package is created for each case
#      ex: zealpy_landau_.so, zealpy_BOT_.so, ...
#----------------------------------------------------------------------
ext_pyreldisp=[]
for mod in pyreldisp_dic:
    pyreldisp_files = zeal_files_list1 + \
        pyreldisp_dic[mod] + zeal_files_list2
    pyreldisp_files.insert(0,'zeal_%s.pyf'%mod)

    ext_pyreldisp.append(Extension(name='zealpy_%s_'%mod,
                           sources = pyreldisp_files,
                           extra_link_args=extra_link_args,
                           include_dirs=['.']))
#end for
			
setup(name        = "pyreldisp_",
      version     = '0.1',
      description = "python interface for ZEAL code",
      author      = "Eric Sonnendrucker",
      author_email= 'sonnen@ipp.mpg.de',
      url         = 'http://www-irma.u-strasbg.fr/~sonnen',
      ext_modules = ext_pyreldisp
      )

