import numpy as np
import os
import sys
from scipy.interpolate import griddata

#(abs((0.5*$3+$4)/(0.5*$5+$6)

def operation(f):
  #return max(f)
  #return max(np.sqrt(f[:,12]))
  return max(abs((0.5*f[:,1]+f[:,2])/(0.5*f[:,3]+f[:,4])-1.))
  
def extension_name():
  return "nrj"  