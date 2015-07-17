#python3 my_test.py
import os
import sys
import unittest
import itertools
from sets import Set

#sys.path.insert(1, '../')
#sys.path.insert(1, '/Users/mehrenbe/prog/f90nml-0.10.2/')
import f90nml
from f90nml.fpy import f90repr, pybool

num_param = len(sys.argv)

if (num_param == 2):
  nml_file = str(sys.argv[1])+'.nml'
  nml_dict = f90nml.read(nml_file)
  param_file = str(sys.argv[1])+'.param'
  if(os.path.isfile(param_file) is False):
    print("file "+str(param_file)+" does not exist")
    sys.exit()
  infile = open(param_file,'r')
  num_loops = 0
  for line in infile:
    tmp = line.split()
    if(len(tmp)!=0):
      if(tmp[0] == '#'):
        num_loops = num_loops+1
  infile.close()
  #print("num_loops="+str(num_loops))

  infile = open(param_file,'r')    
  num_loops_loc = -1
  num_run = [0 for i in range(num_loops)]
  run = [0 for i in range(num_loops)]
  param = [0 for i in range(num_loops)]
  for line in infile:
    tmp = line.split()
    if(len(tmp)!=0):
      if(tmp[0] == '#'):
        num_loops_loc = num_loops_loc+1
        num_run[num_loops_loc] = 0
        param[num_loops_loc] = tmp[1:]
        run[num_loops_loc] = []
      if (num_loops_loc>=0) and (tmp[0] != '#'):
        run[num_loops_loc].append(tmp)
        num_run[num_loops_loc] = num_run[num_loops_loc]+1
  infile.close()       
  #print("num_run=",num_run)  
  #print("param=",param)
  #print("run=",run)  
  num_loops_loc = -1
  param_nml = [0 for i in range(num_loops)]
  for param_loc in param:
    num_loops_loc = num_loops_loc+1
    tmp=[]
    tmp2=[]
    for p in param_loc:
      for i in nml_dict.keys():
        if p in nml_dict[i]:
          tmp2.append(type(nml_dict[i][p])) 
          tmp.append(i)
    param_nml[num_loops_loc] = tmp
  #print("param_nml=",param_nml)

  #check if param and param_nml have same size
  if(len(param)!=len(param_nml)):
    print("param and param_nml have not same size:",len(param),len(param_nml))
    sys.exit()
  if(len(param)==0):
    print("size param=0")
    sys.exit()      
  for i in range(len(param)):
    param_loc = param[i]
    param_nml_loc = param_nml[i]
    if(len(param_loc)==0):
      print("size param_loc=0")
      sys.exit()
    if(len(param_loc)!=len(param_nml_loc)):
      print("param_loc and param_nml_loc have not same size:",len(param_loc),len(param_nml_loc))
      sys.exit()


  i=[0 for j in range(num_loops)]
  s=-1
  for i in itertools.product(*run):
    for j in range(num_loops):
      for k in range(len(param[j])):
        nml_dict[param_nml[j][k]][param[j][k]] = eval(i[j][k])       
    s=s+1
    nml_dict.write(str(sys.argv[1])+'_'+str(s)+'.nml')
elif  (num_param == 3):
  nml_file = str(sys.argv[1])+'.nml'
  nml_dict = f90nml.read(nml_file)
  param_file = str(sys.argv[2])+'.param'
  if(os.path.isfile(param_file) is False):
    print("file "+str(param_file)+" does not exist")
    sys.exit()
  #count the number of loops
  infile = open(param_file,'r')
  num_loops = 0
  for line in infile:
    tmp = line.split()
    if(len(tmp)!=0):
      if(tmp[0] == '#'):
        num_loops = num_loops+1
  infile.close()
  infile = open(param_file,'r')    
  num_loops_loc = -1
  num_run = [0 for i in range(num_loops)]
  run = [0 for i in range(num_loops)]
  param = [0 for i in range(num_loops)]
  for line in infile:
    tmp = line.split()
    if(len(tmp)!=0):
      if(tmp[0] == '#'):
        num_loops_loc = num_loops_loc+1
        num_run[num_loops_loc] = 0
        param[num_loops_loc] = tmp[1:]
        run[num_loops_loc] = []
      if (num_loops_loc>=0) and (tmp[0] != '#'):
        run[num_loops_loc].append(tmp)
        num_run[num_loops_loc] = num_run[num_loops_loc]+1
  infile.close()       
  #print("num_run=",num_run)  
  #print("param=",param)
  #print("run=",run)  
  #sys.exit()
  num_loops_loc = -1
  param_nml = [0 for i in range(num_loops)]
  for param_loc in param:
    num_loops_loc = num_loops_loc+1
    tmp=[]
    tmp2=[]
    for p in param_loc:
      for i in nml_dict.keys():
        if p in nml_dict[i]:
          tmp2.append(type(nml_dict[i][p])) 
          tmp.append(i)
    param_nml[num_loops_loc] = tmp
  #print("param_nml=",param_nml)
  if(len(param)!=len(param_nml)):
    print("param and param_nml have not same size:",len(param),len(param_nml))
    sys.exit()
  if(len(param)==0):
    print("size param=0")
    sys.exit()      
  for i in range(len(param)):
    param_loc = param[i]
    param_nml_loc = param_nml[i]
    if(len(param_loc)==0):
      print("size param_loc=0")
      sys.exit()
    if(len(param_loc)!=len(param_nml_loc)):
      print("param_loc and param_nml_loc have not same size:",len(param_loc),len(param_nml_loc))
      sys.exit()
  i=[0 for j in range(num_loops)]
  s=-1
  setdirname=Set()
  for i in itertools.product(*run):
    filname=str(sys.argv[1])
    dirname=str(sys.argv[2])
    for j in range(num_loops):
      for k in range(len(param[j])):
        nml_dict[param_nml[j][k]][param[j][k]] = eval(i[j][k])
      #if(len(i[j])>len(param[j])):
      #print("size:",j,len(i[j]),len(param[j]))
      if(len(i[j])==len(param[j])+1):      
        #dirname+=str(i[j][len(param[j])])
        #print("val:",i[j][len(param[j])])
        check=0
        check_count=0
        for ell in i[j][len(param[j])]:
          if(ell=='/'):
            check=check+1
            check_val=check_count 
          check_count +=1
        if(check!=1):
          print("i[j][len(param[j])]=",i[j][len(param[j])])
          print("check=",check)
          print("check should be =1")
          sys.exit()
        #print("check_val=",check_val)           
        before=i[j][len(param[j])][0:check_val]
        after=i[j][len(param[j])][check_val+1:]
        dirname+=before
        filname+=after
        #print("before=",before)
        #print("after=",after)
      else:
        print("i[j]=",i[j])
        print("len(i[j])=",len(i[j]))
        print("param[j]=",param[j])
        print("len(param[j])=",len(param[j]))
        print("we should have: len(i[j])==len(param[j])+1")
        sys.exit()      
    #print("dirname=",dirname,s)
    #list_dirname.append(dirname)
    filename=dirname+'/'+filname
    setdirname.add(dirname)
    if not os.path.exists(os.path.dirname(filename)):
      os.makedirs(os.path.dirname(filename))
    s=s+1
    nml_dict.write(filename+'.nml')
  linux_runs='runs=('
  for i in setdirname:
    print(i)
    linux_runs+=' '+i
  linux_runs+=' )'
  print(linux_runs)  
  #print('list=',list_dirname)
else:  
  print("syntax: python /Users/mehrenbe/prog/f90nml-0.10.2/test/my_test.py filename")
  print("we suppose that filename.nml exists")
  print("we suppose that filename.param exists")
  sys.exit()
    
  #print(nml_dict)
#for f in sys.argv[1:]:
#  print(f)
#  nml = f90nml.read(f)
#  print(nml)
#nml = f90nml.read(sys.argv[1])