from zealpy_landau import zealpy_landau
from pylab import *
import numpy as np

rcParams['interactive']=True

for kmode in arange(0.5,0.0,-.1):
    zp = zealpy_landau(kmode=kmode)
    zp.get_zeros(-10,10,-2,10)

    print 'mode:', kmode 
    print 'zero with largest imaginary part:', \
        zp.zeros[argmax(imag(zp.zeros))]

    plot(np.real(zp.zeros),np.imag(zp.zeros),'.',label='k='+str(kmode))
    axis([-3,3,-2.2,.2])
    title('zeros of dispersion relation')
    legend()

if ( (rcParams['backend']=='TkAgg') or
     (rcParams['backend']=='WXAgg') ):
    show()

raw_input("Press enter") 
