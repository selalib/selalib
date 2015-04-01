import numpy as np
import os
import sys
from scipy.interpolate import griddata

def operation(f):
  #return max(f)
  #return max(np.sqrt(f[:,12]))
  return max(abs(np.sqrt(f[:,9])-np.sqrt(f[:,12]))/np.sqrt(f[:,12]))
  
def extension_name():
  return "l2"  