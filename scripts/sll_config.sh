#!/bin/bash
EXTRA_ARGS=$@
hostname=`hostname -s`
echo "$hostname"
case "${hostname}" in 'irma-suse') #irma-suse (Ubuntu 12.10)
        
        cmake \
          -D CMAKE_BUILD_TYPE:STRING=DEBUG \
          $EXTRA_ARGS \
          ${SELALIB_HOME}

;;
'mp-navaro') 
echo "mac book pro intel"
source /opt/intel/composerxe/bin/compilervars.sh intel64
source /opt/intel/mkl/bin/mklvars.sh intel64
export FC=ifort
export OMPI_FC=ifort
export HDF5_HOME="/opt/intel"
export FFTW_HOME="/opt/intel"
;;
'irma-hpc') #irma-hpc (Debian sid)
echo "*** For Intel compilers ***"
source /etc/bash_completion.d/git
export PS1='\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;31m\]\w\[\033[00m\]\e[0;36m$(__git_ps1 " (%s)")\e[m\$ '
source /opt/intel/composerxe/bin/compilervars.sh intel64
source /opt/intel/mkl/bin/mklvars.sh intel64
export PATH=/opt/local/bin:${PATH}
export FC=ifort
export F77=ifort
export CC=icc
export CXX=icpc
export F9X=ifort
export OMPI_FC=ifort
export OMPI_CC=icc
export OMPI_CXX=icpc
export HDF5_ROOT=/opt/local
export FFTW_ROOT=/opt/local
export CFLAGS='-O3 -xHost -ip'
export CXXFLAGS='-O3 -xHost -ip'
export FFLAGS='-O3 -xHost -ip'
export F90FLAGS='-O3 -xHost -ip'
;;
'irma-spare') # irma-spare (Fedora Core 17) with intel compilers
echo "*** For Intel compilers ***"
source /opt/intel/composerxe/bin/compilervars.sh intel64
source /opt/intel/mkl/bin/mklvars.sh intel64
export FC=ifort
export F77=ifort
export CC=icc
export CXX=icpc
export MPICH_FC=ifort
export MPICH_F90=ifort
export MPICH_CC=icc
export MPICH_CXX=icpc
export HDF5_ROOT=/opt/local
export FFTW_ROOT=/opt/local
mpd &
;;
'hpc')
export PATH=~navaro/local/bin:${PATH}
module load compilers/intel11
module load libs/mkl11 
export CC=icc
export FC=ifort
export HDF5_ROOT=~navaro/local/bin
export FFTW_ROOT=~navaro/local
export OMPI_CC=icc
export OMPI_FC=ifort
export OMPI_MCA_btl=self,sm
module list
source /etc/bash_completion.d/git
PS1='[\u@\h \W$(__git_ps1 " (%s)")]\$ '
export GIT_PS1_SHOWDIRTYSTATE=true
export GIT_PS1_SHOWUPSTREAM="auto"
export LD_LIBRARY_PATH=${HOME}/local/lib:\${LD_LIBRARY_PATH}
;;
esac

echo "FC:${FC}"
echo "OMPI_FC:${OMPI_FC}"
echo "HDF5_ROOT:${HDF5_ROOT}"
echo "FFTW_ROOT:${FFTW_ROOT}"
