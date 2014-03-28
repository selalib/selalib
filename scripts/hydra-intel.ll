# @ shell=/bin/bash
# @ error   = selalib.err.$(jobid)
# @ output  = selalib.out.$(jobid)
# @ job_type = parallel
# @ node_usage = not_shared
# @ node = 4
# @ tasks_per_node = 16
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = 0:10:00
# @ notification = complete
# @ notify_user = $(user)@rzg.mpg.de
# @ queue
cd /ptmp/pin/all
module load intel
module load mpi.ibm
module load cmake
module load hdf5-mpi
module load fftw
module load git
make NightlyTest
