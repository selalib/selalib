# @ shell=/bin/bash
#
# @ error   = job.err.$(jobid)
# @ output  = job.out.$(jobid)
# @ job_type = parallel
# @ environment = COPY_ALL
# @ node_usage= not_shared
# @ node = 1
# @ tasks_per_node = 1
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = 00:05:00
# @ notification = complete
# @ queue

#
# run the program
#
OMP_NUM_THREADS=16 
export OMP_NUM_THREADS

poe ./fmurge 1000 3  > prog.out 
#####################################################
