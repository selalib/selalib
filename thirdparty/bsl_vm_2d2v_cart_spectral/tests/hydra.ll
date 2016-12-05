# @ shell=/bin/bash
# @ error   = vm4d.err.$(jobid)
# @ output  = vm4d.out.$(jobid)
# @ job_type = parallel
# @ node_usage = not_shared
# @ node = 4
# @ tasks_per_node = 16
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = 20:00:00
# @ notification = complete
# @ notify_user = $(user)@rzg.mpg.de
# @ queue
module load intel
module load mpi.ibm
module load hdf5-mpi
module load fftw
cd /ptmp/${USER}/vm4d/runs
