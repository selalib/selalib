import os
from pylab import *
import numpy as n

N_x1_tab=arange(32,257,32)#[64,128]
N_x2_tab=N_x1_tab#[32,64]
execu=1


if execu==1:
for N_x1 in N_x1_tab:
#print N_x1
#for N_x2 in N_x2_tab:
strg="echo \"&param N_x1="+str(N_x1)+",N_x2="+str(N_x1)
strg=strg+"/\" >to;./rotation_CSL2D <to >>tf"
print strg
os.system(strg)
os.system("mv tf resb.dat")

