from zealpy_BOT import zealpy_BOT
from pylab import *
import numpy as np

rcParams['interactive']=True

kmode = 0.5

xmin = -2.
xmax = 2.
ymin = 0.0001
ymax = 1.

for kmode in arange(0.5,0.0,-.1):
    zp = zealpy_BOT(
        kmode  = kmode,
        a_bump = 0.4,
        v0     = 2.0,
        sigmab = 0.1,
        kbT    = 1.0,
        l0ld   = 18.0)
    zp.get_zeros(xmin,xmax,ymin,ymax)

    print 'mode:', kmode 
    print 'zero with largest imaginary part:', \
        zp.zeros[argmax(imag(zp.zeros))]

    plot(np.real(zp.zeros),np.imag(zp.zeros),'.',label='k='+str(kmode))
    axis([-1.,1.,-1.,1.])
    title('zeros of dispersion relation')
    legend()
#end for

if ( (rcParams['backend']=='TkAgg') or
     (rcParams['backend']=='WXAgg') ):
    show()

raw_input("Press enter") 
