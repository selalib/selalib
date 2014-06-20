#!/bin/bash
source /home/configfiles/bashrc.default
module purge
module load batch/slurm
module load compilers/intel13
module load mpi/openmpi-1.6.i13.threaded
module load libs/mkl13
module load libs/hdf5-1.8
module list
source /etc/bash_completion.d/git
PS1='\h:\w\e[0;36m$(__git_ps1 " (%s)")\e[m\$ '
export GIT_PS1_SHOWDIRTYSTATE=true
export GIT_PS1_SHOWUPSTREAM="auto"
export HDF5_ROOT=/usr/local/hdf5-1.8.5
export FFTW_ROOT=/usr/local/fftw-3.3.2
source /opt/intel/composer_xe_2013.2.146/mkl/bin/mklvars.sh intel64
cd /workdir/math/navaro/all
make NightlyUpdate
make NightlyConfigure
make NightlyBuild
make NightlyTest
make NightlySubmit
