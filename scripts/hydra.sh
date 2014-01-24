cd /ptmp/pin/all
module purge
module load intel
module load mpi.ibm
module load cmake
module load hdf5-mpi
module load fftw
module load git
make Experimental
make NightlyUpdate
make NightlyConfigure
make NightlyBuild
/usr/bin/llsubmit ${HOME}/selalib/scripts/hydra-intel.ll
make NightlySubmit
