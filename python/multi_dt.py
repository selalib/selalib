#python /Users/mehrenbe/scratch/splitvp/splitvp.py
#cd /Users/mehrenbe/scratch/splitvp/irma/run3
#python /Users/mehrenbe/selalib_mac/selalib/python/multi_nml.py splitvp
#rsync -av splitvp_*.nml mehrenbe@irma-hpc2:/scratch/mehrenbe/splitvp/run5/.
#rsync -av run.sh mehrenbe@irma-hpc2:/scratch/mehrenbe/splitvp/run6/.
#rsync -av -e ssh irmahpc2:/scratch/mehrenbe/splitvp/run7/thdiag_*.dat /Users/mehrenbe/scratch/splitvp/irma/run7/.
#python /Users/mehrenbe/scratch/splitvp/run13/multi_thdiag.py thdiag 0 25
#rsync -av -e ssh helios:/csc/workdir2/mehren/splitvp/runA8/thdiag_*.dat /Users/mehrenbe/scratch/splitvp/helios/runA8/. 
#rsync -av -e ssh splitvp_*.nml helios:/csc/workdir2/mehren/splitvp/runA11/.
#rsync -ave ssh splitvp_*.nml curie:/ccc/scratch/cont003/gen7387/mehrenbm/splitvp/run13/.
#rsync -ave ssh $RUNDIR_local/splitvp/run13/splitvp_*.nml curie:$RUNDIR_curie/splitvp/run13/.
#rsync -av -e ssh curie:/ccc/scratch/cont003/gen7387/mehrenbm/splitvp/run11/thdiag_*.dat /Users/mehrenbe/scratch/splitvp/curie/run11/. 
#python $sll_py/multi_dt.py
import numpy as np
import os
import sys

dict = {}
dict["dt_min"] = 0.1
dict["dt_max"] = 4.
dict["num_dt"] = 100
dict["T_max"] = 60.

for farg in sys.argv[1:]:
  (arg,val) = farg.split("=")
  dict[arg] = val

dict["dt_min"] = float(dict["dt_min"])
dict["dt_max"] = float(dict["dt_max"])
dict["num_dt"] = int(dict["num_dt"])
dict["T_max"] = float(dict["T_max"])



#generation of splitvp.param
A=np.logspace(-np.log(dict["dt_max"])/np.log(2),-np.log(dict["dt_min"])/np.log(2),\
num=dict["num_dt"],base=2)
nbstep=np.floor(dict["T_max"]*A+1)
print("# dt number_iterations")
for i in range(dict["num_dt"]):
  print("%1.20g %d" % (1/A[i],nbstep[i]))



