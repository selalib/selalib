#--------------------------------------------------------------------
# Find and plot the frequencies and the growth rate
#  for all (m,n) modes for the drift-kinetic 4D slab case
#--------------------------------------------------------------------
import numpy as np
from   zealpy_gyrokin_anal import zealpy_gyrokin_anal
import matplotlib.pyplot as mpp
#import pyreldisp_plot as pyp

#------------------------------------
# INPUT DATAS (defined as GYSELA)
#------------------------------------
Zi           = 1
#--> For geometry definition
NNr          = 40
rhomin       = 0.00001
rhomax       = 1.
minor_radius = 14.5
aspect_ratio = 16.5
#--> For analytical density and temperature profiles definition 
kappaTi      = 66.
kappaTe      = 66.
kappan       = 13.2
deltarTi     = 0.1
deltarTe     = 0.1
deltarn      = 0.2
#--> Set of considered modes (m,n) 
#-->  for (mm_min < m < mm_max , nn_min < n < nn_max)
mm_min       = 1
mm_max       = 10 #32
nn_min       = 1
nn_max       = 4
#--> Choice of the box where zeros will be search
xmin         = -1.
xmax         = 1.
ymin         = 0.0001
ymax         = 0.1
iota = 0.001

#------------------------------------------------------------
# Normalization to obtain the same profiles than in GYSELA
#------------------------------------------------------------
#--> Geometry definition
R0   = minor_radius*aspect_ratio
rmin = rhomin*minor_radius
rmax = rhomax*minor_radius
Lr   = np.abs(rmax-rmin)
Lz   = 2.*np.pi*R0
#--> Same normalization than GYSELA for profiles 
deltarTi = deltarTi*Lr
deltarTe = deltarTe*Lr
deltarn  = deltarn*Lr
invLTi   = kappaTi/R0
invLTe   = kappaTe/R0
invLn    = kappan/R0


#------------------------------------------------------------
# Dispersion relation solving
#------------------------------------------------------------
#--> Constructor for using dispersion relation for 
#-->   gyrokinetic with analytic n and T profiles
zp = zealpy_gyrokin_anal(
    Zi=Zi,
    NNr=NNr,rmin=rmin,Lr=Lr,
    invLTi=invLTi,deltarTi=deltarTi,
    invLTe=invLTe,deltarTe=deltarTe,
    invLn0=invLn,deltarn0=deltarn,iota=iota)

nb_mm = mm_max - mm_min + 1
nb_nn = nn_max - nn_min + 1
mostunstable = np.zeros([nb_mm,nb_nn],dtype=complex)
#--> Search zeros of the previous dispersion relation
for nn in np.arange(nn_min,nn_max+1,1):
    indx_nn = nn-nn_min
    for mm in np.arange(mm_min,mm_max+1,1):
        zp.zeros = []
        indx_mm = mm-mm_min
        kpar    = float(nn) * 2*np.pi/Lz
        ktheta  = float(mm)
        print '==> Search zeros of mode (m,n) = (%s,%s)'%(str(mm),str(nn))

        #---> Initialisation of the dispersion relation for (m,n) mode
        zp.init_rel_disp(ktheta=ktheta,kpar=kpar)

        #---> Find all the zeros contained is the box defined before
        zp.get_zeros(xmin,xmax,ymin,ymax)

        if (np.size(zp.zeros)!=0):
            #---> Select the most unstable modes
            mostunstable[indx_mm,indx_nn] = \
                zp.zeros[np.argmax(np.imag(zp.zeros))]
            print '----> Zero with largest imaginary part:', \
                mostunstable[indx_mm,indx_nn]
            print '----> All zeros found:', zp.zeros
            print ' '
        else:
            print '---> BECAREFUL: Zeros not found for (m,n) = (%s,%s)'%(str(mm),str(nn)) 
    #end for
#end for

#------------------------------------------------------------
# Figures
#------------------------------------------------------------
#--> Figure 1: Radial density, temperature and eta profiles
rmesh = zp.params.dic['rmesh']
Ti    = zp.params.dic['Ti']
Te    = zp.params.dic['Te']
n0    = zp.params.dic['n0']
eta   = zp.eta
pyp.plot_nT_profiles(rgrid=rmesh,Ti=Ti,Te=Te,n0=n0,eta=eta)


#--> Figure 2: Plot growth rate versus m for differente values of n
pyp.plot_gammam(mm_min,mm_max,nn_min,nn_max,mostunstable,ifig=2)

#---> Display figures at screen
if ( (mpp.rcParams['backend']=='TkAgg') or
     (mpp.rcParams['backend']=='WXAgg') ):
    mpp.show()

raw_input("Press enter")

