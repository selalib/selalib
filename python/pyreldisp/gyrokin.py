import numpy as np
from zealpy_gyrokin_anal import zealpy_gyrokin_anal
from pylab import *

# Gysela-like parameters
Zi   = 1
NNr  = 100

rhomin       = 0.00001
rhomax       = 1.
minor_radius = 14.5
aspect_ratio = 16.5

R0           = minor_radius*aspect_ratio
rmin         = rhomin*minor_radius
rmax         = rhomax*minor_radius
Lr           = np.abs(rmax-rmin)
Lz           = 2.*pi*R0

kappaTi      = 66.
kappaTe      = 66.
kappan       = 13.2
deltarTi     = 0.1*Lr
deltarTe     = 0.1*Lr
deltarn      = 0.2*Lr

# -- after this line, computation of the init profiles --
invLTi = kappaTi/R0
invLTe = kappaTe/R0
invLn  = kappan/R0
zp     = zealpy_gyrokin_anal(Zi=Zi,
                    NNr=NNr,rmin=rmin,Lr=Lr,
                    invLTi=invLTi,deltarTi=deltarTi,
                    invLTe=invLTe,deltarTe=deltarTe,
                    invLn0=invLn,deltarn0=deltarn)

# mode considere
mm     = 5
nn     = 1
kpar   = float(nn) * 2*np.pi/Lz
ktheta = float(mm)
print 'mode',mm,nn
zp.init_rel_disp(ktheta=ktheta,kpar=kpar)

# definition de la boite dans laquelle on cherche les zeros
xmin = -1.
xmax = 1.
ymin = 0.0001
ymax = 0.1
zp.get_zeros(xmin,xmax,ymin,ymax)

mostunstable= zp.zeros[argmax(imag(zp.zeros))]
print 'zero with largest imaginary part:', mostunstable
print 'zeros', zp.zeros
figure(1)

plot(np.real(zp.zeros),np.imag(zp.zeros),'.',
     label='m='+str(mm)+' n='+str(nn))
axis([xmin,xmax,ymin,ymax])
title('zeros of dispersion relation')
legend()
draw()

# Calcul du vecteur propre (phi) correspondant au mode le plus instable
zp.function_input_module.kernel(mostunstable)
# rescale eigenvector to get phi from phi
rmesh = zp.params.dic['rmesh']
n0    = zp.params.dic['n0']
phi_eigvect = zp.function_input_module.vecteur[:,0] / \
    np.sqrt(rmesh*n0)
figure(2)
plot(rmesh,abs(phi_eigvect))
title('Profil du mode le plus instable (Phi)')

if ( (rcParams['backend']=='TkAgg') or
     (rcParams['backend']=='WXAgg') ):
    show()

raw_input("Press enter")

