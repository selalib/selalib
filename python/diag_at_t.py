import numpy as np
import os
import sys
from scipy.interpolate import griddata

def operation(f,tval):
  return griddata(f[:,0], f[:,:], [tval], method='cubic')
  
def extension_name():
  return "newval"  