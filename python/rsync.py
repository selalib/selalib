"""
will dot the following transfers:

  synchronisation from $RUNDIR_local/curie/splitvp/run1/splitvp_*.nml 
  to $RUNDIR_curie/splitvp/run1/.

  synchronisation from $RUNDIR_curie/splitvp/run1/thdiag_*.dat 
  to $RUNDIR_local/curie/splitvp/run1/.

note that a lign that begin with "## " is a commented


following type of lines should be in ~/.bash_profile or ~/.bashrc on local computer:

export RUNDIR_curie=/ccc/scratch/cont003/gen7387/mehrenbm
export RUNDIR_local=$HOME/scratch
export sll_py=$HOME/selalib_mac/selalib/python

following line should be in ~/.bash_profile or ~/.bashrc on dist computer:

export sll_py=$HOME/selalib_mac/selalib/python

do not forget to do (if necessary)

source $HOME/.bash_profile


Example of file (here LOCAL_FILES are commented)
Begin of file splitvp_curie.backup (same format should be used)

## here commands of different scripts are recalled
##
## PREPARATION STEP: generate splitvp.nml and splitvp.param in run directories
##
## suppose that run directories are: run8g run10g run11g run12g run13g run14g
## splitvp.nml can be copied from run8f run10f run11f run12f run13f run14f
## splitvp.param should be updated to new test, for example, change initial condition
## for generating several dt use 
## python $sll_py/multi_dt.py dt_min=0.1 dt_max=4 num_dt=100 T_max=60
## adapt job_file.sh and put together with splitvp.param to run8g run10g run11g run12g run13g run14g
## here are the possible commands, once job*.sh splitvp.param 
## have been updated on $RUNDIR_local/splitvp/curie
##
## old_runs=( run8f run10f run11f run12f run13f run14f )
## runs=( run8g run10g run11g run12g run13g run14g )
## for (( i=0; i<${#runs[@]};i++ ));do cp ${old_runs[$i]}/splitvp.nml ${runs[$i]}/.;done 
## for i in ${runs[@]};do cp job*.sh splitvp.param $i/.;done
##
## FIRST STEP: generate the splitvp_num.nml in each directory
##  
## python $sll_py/map2dir.py "python $sll_py/multi_nml.py splitvp" run8g run10g run11g run12g run13g run14g
##
## SECOND STEP: copy to dist computer (here curie; uncomment LOCAL_FILES; comment DIST_FILES)
## 
## python $sll_py/rsync.py $RUNDIR_local/splitvp/curie/splitvp_curie
##
## THIRD STEP: do the submission on dist computer
##
## check before that ##
## python $sll_py/map2dir.py "ccc_msub job*.sh" run8g run10g run11g run12g run13g run14g
##
## FOURTH STEP: copy from dist computer (here curie; uncomment LOCAL_FILES; comment DIST_FILES)
##
## python $sll_py/rsync.py $RUNDIR_local/splitvp_curie
##
## FIFTH STEP: do a post-processing diag
##
## python $sll_py/map2dir.py "python $sll_py/multi_thdiag.py thdiag 0 99 6 15 30 50" run8g run10g run11g run12g run13g run14g 
##

# DIST_ROOT
curie

# LOCAL_ROOT
local

# RUN_NAME
splitvp


# RUN_DIRS
run8g
run10g
run11g
run12g
run13g
run14g


# LOCAL_FILES
## splitvp_*.nml
## splitvp.param
## job*.sh


# DIST_FILES
thdiag_*.dat
End of file splitvp_.backup

"""

import os
import sys

num_param = len(sys.argv)

#print("num_param=",num_param)

if(num_param != 2):
  print("Usage python rsync.py filename")
  print("filename.backup should be present")
  sys.exit()

backup_file = sys.argv[1]+".backup"

if not os.path.isfile(backup_file):
  print(backup_file+" does not exist")
  sys.exit()

file_id = open(backup_file,'r')
num_loops = 0
for line in file_id:
  tmp = line.split()
  if(len(tmp)!=0):
    if(tmp[0] == '#'):
      num_loops = num_loops+1
file_id.close()

#print("num_loops="+str(num_loops))

file_id = open(backup_file,'r')    
num_loops_loc = -1
num_run = [0 for i in range(num_loops)]
run = [0 for i in range(num_loops)]
param = [0 for i in range(num_loops)]
for line in file_id:
  tmp = line.split()
  if(len(tmp)!=0):
    if(tmp[0] == '#'):
      num_loops_loc = num_loops_loc+1
      num_run[num_loops_loc] = 0
      param[num_loops_loc] = tmp[1:]
      run[num_loops_loc] = []
    if (num_loops_loc>=0) and (tmp[0] != '#') and (tmp[0] != '##'):
      run[num_loops_loc].append(tmp)
      num_run[num_loops_loc] = num_run[num_loops_loc]+1
file_id.close()       
#print("num_run=",num_run)  
#print("backup_keys=",param)
#print("run=",run)

backup_dict = dict()
for i in range(num_loops):
  backup_dict[param[i][0]] = [run[i][j][0] for j in range(num_run[i])]

#print("backup_dict=",backup_dict)

set_backup_keys = {"DIST_ROOT","LOCAL_ROOT","RUN_NAME","RUN_DIRS","LOCAL_FILES","DIST_FILES"}

if( not "DIST_ROOT" in backup_dict.keys()):
  print("DIST_ROOT is not present")
  sys.exit()

if( not "LOCAL_ROOT" in backup_dict.keys()):
  print("LOCAL_ROOT is not present")
  sys.exit()

if( not "RUN_NAME" in backup_dict.keys()):
  print("RUN_NAME is not present")
  sys.exit()

if( not "RUN_DIRS" in backup_dict.keys()):
  print("RUN_DIRS is not present")
  sys.exit()

if( not "LOCAL_FILES" in backup_dict.keys()):
  print("LOCAL_FILES is not present")
  sys.exit()

if( not "DIST_FILES" in backup_dict.keys()):
  print("DIST_FILES is not present")
  sys.exit()


if set(backup_dict.keys()) != set_backup_keys:
  print("backup_keys mismatch")
  print(set(backup_dict.keys()))
  print(set_backup_keys)
  sys.exit()

#num_local_files = len(backup_dict["LOCAL_FILES"])  
#num_dist_files = len(backup_dict["DIST_FILES"])


#print(len(backup_dict["LOCAL_FILES"]))

#get environment variables

rundir_local=os.environ.get("RUNDIR_"+backup_dict["LOCAL_ROOT"][0])
if(rundir_local==None):
  print("please define RUNDIR_"+backup_dict["LOCAL_ROOT"][0]+" as env variable")
  sys.exit()
rundir_dist=os.environ.get("RUNDIR_"+backup_dict["DIST_ROOT"][0])
if(rundir_dist==None):
  print("please define RUNDIR_"+backup_dict["DIST_ROOT"][0]+" as env variable")
  sys.exit()

if(len(backup_dict["LOCAL_FILES"])>=1):
  str="rsync -avRe ssh"
  for j in backup_dict["RUN_DIRS"]:
    for i in backup_dict["LOCAL_FILES"]:
      str=str+" "+rundir_local+"/"\
      +backup_dict["RUN_NAME"][0]+"/"+backup_dict["DIST_ROOT"][0]+"/./"\
      +j+"/"\
      +i
    

  str = str+" "+backup_dict["DIST_ROOT"][0]+":"+rundir_dist+"/"\
  +backup_dict["RUN_NAME"][0]+"/"

  #print(str)
  #os.system("echo "+str)
  os.system(str)


if(len(backup_dict["DIST_FILES"])>=1):
  str="rsync -avRe ssh "+backup_dict["DIST_ROOT"][0]+":"+"\'"
  tmp = 0
  for j in backup_dict["RUN_DIRS"]:
    for i in backup_dict["DIST_FILES"]:
      if(tmp!=0):
        str=str+" "
      str=str+rundir_dist+"/"\
      +backup_dict["RUN_NAME"][0]+"/./"+j+"/"+i
      tmp =tmp+1
  
  str = str+"\' "+rundir_local+"/"\
  +backup_dict["RUN_NAME"][0]+"/"+backup_dict["DIST_ROOT"][0]+"/"    
  
  #print(str)
  #os.system("echo "+str)
  os.system(str)


  

  