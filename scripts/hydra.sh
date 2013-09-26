# @ shell=/bin/bash
# @ error   = selalib.err.$(jobid)
# @ output  = selalib.out.$(jobid)
# @ job_type = parallel
# @ node_usage= not_shared
# @ node = 8
# @ tasks_per_node = 16
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = 1:00:00
# @ notification = complete
# @ notify_user = $(user)@rzg.mpg.de
# @ queue
module purge
module load intel
module load mkl
module load mpi.ibm
module load cmake
module load hdf5-mpi
module load fftw
module load git
export FC=ifort
export CC=icc
export CXX=icpc
cd /ptmp/pin
mkdir selalib-intel
cd selalib-intel; {
cmake ${HOME}/selalib -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=ifort
make NightlyUpdate
make NightlyConfigure
make NightlyBuild
make NightlyTest
make NightlySubmit
}; cd -
rm -rf /ptmp/pin/selalib-intel
