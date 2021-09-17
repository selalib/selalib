#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 13:43:29 2019

@author: Katharina Kormann
"""

import f90nml as f90nml
import os
import subprocess

folder = 'run_directory'
executable = 'sim_pic_vm_1d2v_cart'
sample_file = 'sample.nml'

file_path = folder + sample_file
os.chdir(folder)

# MPI command 
# TODO: MPI command currently not working
cmd = ''#mpirun -np 1 '

# Read in the sample file
nml = f90nml.read(file_path)

# Change values of the nml structure
# Vector defining the various parameter for which a simulation should be run
sdeg_vec = [1, 2, 3, 4 ]
dfem_vec = [-1,1,2,3,4,5]
propagators = ['disgrad', 'symplectic']

for sdeg in sdeg_vec :
    for dfem in dfem_vec :
        for propagator in propagators :
            # Change the values of the nml file
            nml['pic_params']['spline_degree'] = sdeg
            nml['grid_dims']['degree_fem'] = dfem
            nml['pic_params']['splitting_case'] = 'splitting_'+propagator
            prefix = propagator +'_f'+str(dfem)+'_s' + str(sdeg)
            nml['output']['file_prefix'] = prefix
            file_out = folder + prefix + '.nml'
            # Remove old file (if existing)
            if os.path.exists(file_out):
                os.remove(file_out)
            # Write the nml file
            f90nml.write(nml, file_out)
            # Finally submit the simulation
            print([executable, file_out])
            out = subprocess.Popen([executable, file_out], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout,stderr = out.communicate()
            print(stdout)
            print (stderr)