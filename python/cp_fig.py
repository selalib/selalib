"""
will do the following print:

  cp SOURCE_ROOT/src_file1 DIST_ROOT/dist_file1 
  cp SOURCE_ROOT/src_file2 DIST_ROOT/dist_file2  
  ...
  cp SOURCE_ROOT/src_fileN DIST_ROOT/dist_fileN  
  
where following examples of lines are stored in sys.argv[1].backup

# SOURCE_ROOT
/Users/mehrenbe/scratch/curvisl/cg/cemracs14/irmahpc2/

# DIST_ROOT
/Users/mehrenbe/selalib_mac/sll_docs/curvisl/cemracs2014/Figures/ 

# FILES
src_file1 dist_file1
src_file2 dist_file2
...
src_fileN dist_fileN

note that a lign that begin with "## " is a commented
in the file argv[1].backup

"""

import os
import sys

num_param = len(sys.argv)

#print("num_param=",num_param)

if(num_param != 2):
  print("Usage python cp_fig.py filename")
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

check_backup_keys = [["SOURCE_ROOT"],["DIST_ROOT"],["SOURCE_FILES","DIST_FILES"]]

if(check_backup_keys != param):
  print("bad sys.argv[1] file")
  sys.exit()

if(len(num_run) !=3 ):
  print("bad size for num_run")
  sys.exit()

if(num_run[0]!=1):
  print("bad size for num_run[0]")
  sys.exit()
if(num_run[1]!=1):
  print("bad size for num_run[1]")
  sys.exit()
if(num_run[2]<1):
  print("bad size for num_run[2]")
  sys.exit()

if(len(run) !=3 ):
  print("bad size for run")
  sys.exit()

if(len(run[0])!=1):
  print("bad size for run[0]")
  sys.exit()
if(len(run[0][0])!=1):
  print("bad size for run[0][0]")
  sys.exit()

if(len(run[1])!=1):
  print("bad size for run[1]")
  sys.exit()
if(len(run[1][0])!=1):
  print("bad size for run[1][0]")
  sys.exit()

if(len(run[2])<1):
  print("bad size for run[2]")
  sys.exit()


for val in run[2]:
  if(len(val) !=2):
    print("bad size for val in run[2]")
    sys.exit()


source_root = run[0][0][0]
dist_root = run[1][0][0]

#print(source_root)
#print(dist_root)


for val in run[2]:
  print "cp "+source_root+val[0]+" "+dist_root+val[1]

sys.exit()


backup_dict = dict()
for i in range(num_loops):
  backup_dict[param[i]] = [run[i][j] for j in range(num_run[i])]

#print("backup_dict=",backup_dict)

set_backup_keys = {["SOURCE_ROOT"],["DIST_ROOT"],["SOURCE_FILES","DIST_FILES"]}

if( not ["SOURCE_ROOT"] in backup_dict.keys()):
  print("SOURCE_ROOT is not present")
  sys.exit()


if( not ["DIST_ROOT"] in backup_dict.keys()):
  print("DIST_ROOT is not present")
  sys.exit()

if( not ["SOURCE_FILES","DIST_FILES"] in backup_dict.keys()):
  print("SOURCE_FILES is not present")
  sys.exit()





if set(backup_dict.keys()) != set_backup_keys:
  print("backup_keys mismatch")
  print(set(backup_dict.keys()))
  print(set_backup_keys)
  sys.exit()

print(backup_dict)

sys.exit()

#num_local_files = len(backup_dict["LOCAL_FILES"])  
#num_dist_files = len(backup_dict["DIST_FILES"])


#print(len(backup_dict["LOCAL_FILES"]))



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


