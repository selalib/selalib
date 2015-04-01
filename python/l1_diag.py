import numpy as np
import os
import sys
from scipy.interpolate import griddata

#(abs((0.5*$3+$4)/(0.5*$5+$6)

def operation(f):
  #return max(f)
  #return max(np.sqrt(f[:,12]))
  return max(abs(f[:,8]-f[:,11])/f[:,11])
  
def extension_name():
  return "l1"  