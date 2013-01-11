import numpy as np
import matplotlib.pylab as plb
import matplotlib.pyplot as plt
import utils_func as ut 
import h5py
import glob

#*** Mesh loading ***
fx1 = h5py.File("polar_mesh-x1.h5","r")
x1  = fx1['x1']
nx1 = np.size(x1)
fx2 = h5py.File("polar_mesh-x2.h5","r")
x2  = fx2['x2']
nx2 = np.size(x2)

#*** f2D initial time loading ***
ff2D_t0 = h5py.File('f0001-values.h5')
f2D_t0  = np.array(ff2D_t0['values'])

#*** look at all the existing result files ***
fnames_data = glob.glob('f*-values.h5')
fnames_data.sort()

#***
nbfiles     = np.size(fnames_data)
nbstep_diag = 1   #10
ifil_beg    = 50  #0
ifil_end    = nbfiles+nbstep_diag-1

#*** computation of the min and max values ***
min_f2D = 10000.
max_f2D = -10000.
for ifil in range(ifil_beg,ifil_end,nbstep_diag):
    #--> f2D loading
    filename = fnames_data[ifil]
    ff2D     = h5py.File(filename)
    f2D      = np.array(ff2D['values'])
    min_f2D  = min(min_f2D,f2D.min())
    max_f2D  = max(max_f2D,f2D.max())
#end for
print 'min f2D = '+str(min_f2D)
print 'max f2D = '+str(max_f2D)
fig = plt.figure()
for ifil in range(ifil_beg,ifil_end,nbstep_diag):
    #--> f2D loading
    filename = fnames_data[ifil]
    print 'filename = ' + filename
    ff2D = h5py.File(filename)
    f2D  = np.array(ff2D['values'])
    #--> f2D plot
    if ( ifil == ifil_beg ):
        fig.clf()
    else:
        ax1.cla();
    ax1 = fig.add_subplot(111)
    ax1.axis('equal')
    p1  = ax1.contourf(x1,x2,f2D,100)
    p1.set_clim((min_f2D,max_f2D));
    plt.title(filename)
    ax1.axis('off')
    if (ifil == ifil_beg):
        fig.colorbar(p1);
    fig.canvas.draw();
#end for

