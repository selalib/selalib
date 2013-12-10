import os
from pylab import *
import numpy as n

N_x1_tab=[31] #arange(32,65,32)#[64,128]
N_x2_tab=N_x1_tab#[32,64]
execu=1


if execu==1:
  for N_x1 in N_x1_tab:
    #print N_x1
    #for N_x2 in N_x2_tab:
    strg="echo \"&param N_eta1="+str(N_x1)+",N_eta2="+str(N_x1)
    strg=strg+"/\" >to;./bin/test_cg_curvilinear <to >>tf"
    print strg
    os.system(strg)

