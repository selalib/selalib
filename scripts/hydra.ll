# @ shell=/bin/bash
# @ error   = selalib.err.$(jobid)
# @ output  = selalib.out.$(jobid)
# @ job_type = parallel
# @ node_usage = shared
# @ node = 8
# @ tasks_per_node = 8
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = 0:10:00
# @ notification = complete
# @ notify_user = $(user)@rzg.mpg.de
# @ queue
cd /ptmp/pin/all
make Experimental
