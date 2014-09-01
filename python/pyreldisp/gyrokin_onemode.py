import numpy as np
from zealpy_gyrokin_anal import zealpy_gyrokin_anal
import matplotlib.pyplot as mpp
from pylab import *

#------------------------------------
# INPUT DATAS (defined as GYSELA)
#------------------------------------
Zi           = 1
#--> For geometry definition
NNr          = 200
rhomin       = 0.006896551724 #0.00001
rhomax       = 1.
minor_radius = 14.5
aspect_ratio = 16.53849335 #16.5
#--> For analytical density and temperature profiles definition 
kappaTi      = 66.
kappaTe      = 66.
kappan       = 13.2
deltarTi     = 0.1
deltarTe     = 0.1
deltarn      = 0.2
#--> Choice of the box where zeros will be search
xmin         = -1.
xmax         = 1.
ymin         = 0.0001
ymax         = 0.01 #0.1
iota = 0.
B0 = -1.
#------------------------------------------------------------
# Normalization to obtain the same profiles than in GYSELA
#------------------------------------------------------------
#--> Geometry definition
R0   = minor_radius*aspect_ratio
rmin = rhomin*minor_radius
rmax = rhomax*minor_radius
Lr   = np.abs(rmax-rmin)
Lz   = 2.*pi*R0
#--> Same normalization than GYSELA for profiles 
deltarTi = deltarTi*Lr
deltarTe = deltarTe*Lr
deltarn  = deltarn*Lr
invLTi   = kappaTi/R0
invLTe   = kappaTe/R0
invLn    = kappan/R0
iota     = iota/R0

#------------------------------------------------------------
# Dispersion relation solving
#------------------------------------------------------------
#--> Constructor for using dispersion relation for 
#-->   gyrokinetic with analytic n and T profiles
zp     = zealpy_gyrokin_anal(Zi=Zi,
                    NNr=NNr,rmin=rmin,Lr=Lr,
                    invLTi=invLTi,deltarTi=deltarTi,
                    invLTe=invLTe,deltarTe=deltarTe,
                    invLn0=invLn,deltarn0=deltarn,iota=iota,B0=B0)

#--> Choice of the mode (m,n)
print 'B0=',B0
mm_choice = raw_input('Poloidal mode m ?: ')
nn_choice = raw_input('Toroidal mode n ?: ')
ifig = 1
while ( (mm_choice != '') | (nn_choice != '') ) :
    zp.zeros = []
    mm  = int(mm_choice)
    nn  = int(nn_choice)

    print '==> Search zeros of mode (m,n) = (%s,%s)'%(str(mm),str(nn))
    kpar    = float(nn) * 2*np.pi/Lz
    ktheta  = float(mm)

    #---> Initialisation of the dispersion relation for (m,n) mode
    zp.init_rel_disp(ktheta=ktheta,kpar=kpar)

    #---> Find all the zeros contained is the box defined before
    zp.get_zeros(xmin,xmax,ymin,ymax)

    #---> Select the most unstable modes
    if (size(zp.zeros)!=0):
        imag_zeros   = np.imag(zp.zeros)
        mostunstable = zp.zeros[np.argmax(imag_zeros)]
        print '----> Zero with largest imaginary part: ', mostunstable
        print '----> All zeros found: ', zp.zeros
        print ' '

        #---> Computation of the eigenvector (phi) corresponding to the
        #--->  most unstable mode
        zp.function_input_module.kernel(mostunstable)
        #-----> Rescaling of the eigenvector to get phi from phi
        rmesh = zp.params.dic['rmesh']
        n0    = zp.params.dic['n0']
        phi_eigvect = zp.function_input_module.vecteur[:,0] / \
            np.sqrt(rmesh*n0)

        #--> Figure 1: zeros of dispersion relation
        figure(ifig)
        plot(np.real(zp.zeros),np.imag(zp.zeros),'.',
             label='m='+str(mm)+' n='+str(nn))
        axis([xmin,xmax,ymin,ymax])
        title('zeros of dispersion relation')
        legend()
        draw()
        ifig = ifig + 1

        #--> Figure 2: Plot radial profile of the most unstable mode
        figure(ifig)
        plot(rmesh,abs(phi_eigvect))
        title('Profile of the most unstable (Phi)')
        draw()
        ifig = ifig + 1

        #---> Display figures at screen
        if ( (mpp.rcParams['backend']=='TkAgg') or
             (mpp.rcParams['backend']=='WXAgg') ):
            mpp.show()
    else:
        print '---> BECAREFUL: Zeros not found for (m,n) = (%s,%s)'%(str(mm),str(nn)) 

    #--> Choice of the mode (m,n)
    mm_choice = raw_input('Poloidal mode m ?: ')
    nn_choice = raw_input('Toroidal mode n ?: ')
#end while

raw_input("Press enter")

