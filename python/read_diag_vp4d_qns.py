import utils_func as ut
import matplotlib as mp
import matplotlib.pyplot as mpp

filename = '/Users/back/selalib/prototype/build/vp4D_diag_d00000.h5'
f_diag = ut.loadHDF5(filename)
print "HDF5 file read = " + filename


fig = mpp.figure(figsize=(18,9))


