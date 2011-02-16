******************************************************
*            GYSELA5D general README file            *
******************************************************

This compressed archive contains all the files listed :
  README         : This file.
  src_public/    : This directory contains all the source files.
  wk_public/     : This working directory contains examples of input data files
  script_public/ : This directory contains some scripts for the 
                   exploiting of GYSELA results

---------------------------------------------------------------------------
Contents
---------------------------------------------------------------------------

* Code Description
------------------
 A. Licence information
 B. General description
 C. Coding and Parallelization

* Building the Code
-------------------
 A. Preliminaries
 B. Compiling

* Running the Code
------------------
 A. A set of 4 SIMPLE test cases
 B. Strong scaling on 512, 1024, 2048, 4096 processors 
 C. Weak scaling on 512, 1024, 2048, 4096 processors 

* Exploiting Results
--------------------


---------------------------------------------------------------------------
* Code Description
------------------
 A. License information

   The code is under CECILL-B licence:
   Copyright Status
   !**************************************************************
   !  Copyright Euratom-CEA
   !  Authors : 
   !     Virginie Grandgirard (virginie.grandgirard@cea.fr)
   !     Chantal Passeron (chantal.passeron@cea.fr)
   !     Guillaume Latu (guillaume.latu@cea.fr)
   !     Xavier Garbet (xavier.garbet@cea.fr)
   !     Philippe Ghendrih (philippe.ghendrih@cea.fr)
   !     Yanick Sarazin (yanick.sarazin@cea.fr)
   !  
   !  This code GYSELA (for GYrokinetic SEmi-LAgrangian) 
   !  is a 5D gyrokinetic global full-f code for simulating 
   !  the plasma turbulence in a tokamak.
   !  
   !  This software is governed by the CeCILL-B license 
   !  under French law and abiding by the rules of distribution 
   !  of free software.  You can  use, modify and redistribute 
   !  the software under the terms of the CeCILL-B license as 
   !  circulated by CEA, CNRS and INRIA at the following URL
   !  "http://www.cecill.info". 
   !**************************************************************

 B General description

   GYSELA5D code is based on a Semi-Lagrangian scheme and solves 5D gyrokinetic ion turbulence in tokamak plasmas.

 C Coding and Parallelization

   This version of gysela5D is implemented in Fortran 90, with some calls in C.
   It is an hybrid code with two levels of parallelism: MPI and OpenMP.
   The main output files are in HDF5 format, so the library HDF5 is needed.

* Building The code
-------------------
A. Preliminaries

   At first, a 'ARCH' environment variable must be initialized with the machine name
   For example in the '.bashrc' file, you should add something like 'export ARCH=hpcff'.

   You need ‘GNU make’ tool to generate the executable.
   To adapt the Makefile for a new platform, modify the make.include (the 'Makefile'
   uses directly this 'make.include' file)
     - add a new part according to this ARCH variable
     - check the Fortran and C compiler, the preprocessor directives, 
     - check for library locations:
        - MPI library
        - HDF5 library
     - check all the compiler flags for appropriate values on your system

   For example, for the HPC-FF machine, here is the part of the make.include file:
        ifeq ($(ARCH), hpcff)
          # --> General path for HDF5 library
          BASEHDF5=/root of HDF5 library
          MAKE = make
          # --> MPI Fortran90 compiler
          F90  = mpif90
          # --> C compiler
          CC  = cc
          # **** FORTRAN COMPILING OPTIONS ****
          ifeq ($(MAKECMDGOALS),debug)
             # --> Debug mode + preprocessor
             F90FLAGS = -fpp -g -check bounds
          else
             # --> preprocessor + optimisation + using of the 'ctime.c' file for CPU time computation
             F90FLAGS = -fpp -O3 -DCTIMER 
          endif
          # **** PREPROCESSOR DIRECTIVE ****
          # --> If -DTIMER is specified more informations of the CPU time are detailed in the output file
          TIMERFLAG = -DTIMER
          # --> If MPI2 is available on your machine, use the directive -DMPI2 
          MPI2FLAG = -DMPI2
          # --> The code is available to run without OpenMP (for this you need to disable the following
                OMPFLAG)
          # --> But you have to keep in mind that GYSELA is parallely optimized with an hybrid
                parallelisation MPI/OpenMP 
          OMPFLAG = -openmp
          # --> C directive
          CFLAGS = -DNORMAL -I/usr/local/include
          # --> HDF5 library paths
          HDF5INCLUDE = -I$(BASEHDF5)/lib
          HDF5LIB     = $(BASEHDF5)/lib/libhdf5_fortran.a $(BASEHDF5)/lib/libhdf5.a
          NOHDF5FLAG  = -DNOHDF5
      endif

B. Compiling 

   You have to go to gysela5D/src_public/ directory.

   The make options are the following:
   make           : by default for generating the executable file gysela.exe'
   make debug     : for compiling with debugging options
   make clean     : for removing all object and module files.
   make distclean : for deleting all : object, module and executable filess.
   make noHDF5    : for compiling without HDF5. In this case there is no HDF5 output files.

   For a standard compilation just type 'make'. If you encounter some difficulties with
   HDF5 library, try 'make clean; make noHDF5'.

* Running the code
------------------
   Data input files are available in gysela5D/wk_public directory
   
   > qsub gysela_<data_name>.job
   This create automatically a directory called RESULTS_<data_name>

   If the execution is OK, the directory RESULTS_<data_name> contains:
   --> The executable 'gysela.exe': a copy of the executable
       generated with the source compiling (in src_public) 

   --> ASCII files: 
         - gysela_res.out : Output file (to follow the evolution of the simulation)
         - gysela_log.out : Log file 
         - gysela_res.err : Error output file (which is empty if no error) 
         - gysela_data_r<restart_num>.out : copy of the input data read by the code
         - gysela_mem.out : Diagnostic on the memory required per node +
                            details of the memory required for each array allocation
         - gysela_CL.out  : Saving of the evolution of quantities concerning the 
	                    conservation laws, as number of ions and electrons,
                            L2norm, entropy, etc ...
         - gysela_rprof.out : Saving of the evolution of certain quantities at 
                              the radial position of the biggest gradients, as 
                              the (0,0) mode of the potential and the flux, etc ...
         - outofdomain_mu<mu_num>.txt : list the iterations where particles go out of 
                                        domain in r or theta directions
         - file_list.out : containing only the number of restart already performed
	 - num_restart.out : the index of the restart file (see futher)


   --> HDF5 files containing information corresponding 
       to the initial state of the simulation:
         - coord_system.h5 : Saving of the coordinates of the system, jacobian, etc...
         - init_state_r<restart_num>.h5 : Saving of the initial state, i.e :
                                          - density and temperature profiles,
                                          - magnetic field, etc ...

   --> directories containing HDF5 result files:
         - conservation_laws/ : OD data concerning the conservation laws
         - rprof/             : 1D radial profiles
         - Phi2D/             : 2D cross-section of the electric potential
         - f2D/               : 2D cross-section of the distribution function       
         - moment3D/          : 3D saving of different moments of the distribution function
	                        (integration in the velocity space)
         - Phi3D/             : 3D saving of the electric potential 
         - f5D/               : 5D saving of the distribution function

    --> The restart files containing the distribution function at time tn and tn-1:
         - NbMPIprocess x gysela.rst.n0.fnm1.p<proc_num>.h5    (for time tn-1)
         - NbMPIprocess x gysela.rst.n0.fn.p<proc_num>.h5      (for time tn)
         - NbMPIprocess x gysela.rst.n1.fnm1.p<proc_num>.h5    (for time tn-1)
         - NbMPIprocess x gysela.rst.n1.fn.p<proc_num>.h5      (for time tn)
        Rk : In the case of predictor-corrector algorithm (leap-frog=.false. in the data file),
             only the time tn is saved

	In both cases, there are two sequences, one indexed by 'n0' and 
	the other one indexed by 'n1'. 
	This is done to avoid problem of read-write, i.e :
         - if the restart files which are read correspond to the sequence 'n0' then 
	   the writing at the end of the simulation is performed on the sequence 'n1' 
	   and the opposite otherwise.
 
* Exploiting the results
------------------------   
Several scripts for exploiting the GYSELA results are in the directory 'script_public'.
 --> GnuPlot scripts (for exploiting ASCII output files): 
      - GPgysela_CL     : read 'gysela_CL.out'
                          + plot the time evolution of the number of ions and
                            electrons, the L2norm and the entropy.
      - GPgysela_energy : read 'gysela_CL.out'
                          + plot the time evolution of the potential and
                            kinetic energies
      - GPgysela_rprof  : read 'gysela_rprof.out'
                          + plot the evolution of Phi(0,0), Phi^2 and of 
                            the flux at a radial position (corresponding
                            to the largest region of gradients)

 --> GnuPlot scripts (for exploiting HDF5 output files): 
      - GPgysela_init   : read 'init_state_r000.h5' 
                          + plot initial temperature and density 
                            radial profiles.
      - GPplotHDF5_var1D : plot any 1D variable contained in an HDF5 file
           example : GPplotHDF5_var1D init_state_r000.h5 Ti 

      - GPplotHDF5_var2D : plot any 2D variable contained in an HDF5 file
           example : GPplotHDF5_var2D Phi2D/Phi2D_d00000.h5 Phirphi

      - GPplotHDF5_rth2D : plot any 2D variable contained in an HDF5 file
	                   on a poloidal cross-section. (coherent only
                           for (r,theta) cross-section saving)
           example : GPplotHDF5_rth2D Phi2D/Phi2D_d00000.h5 Phirth
 
