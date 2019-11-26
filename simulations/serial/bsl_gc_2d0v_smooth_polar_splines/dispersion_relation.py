# See Ronald C Davidson, PHYSICS OF NONNEUTRAL PLASMAS (Chapter 6)
# and section 6.3 of https://doi.org/10.1016/j.jcp.2019.108889
import numpy as np

# doi:10.1140/epjd/e2014-50180-9
#l  = 2
#rm = 4.
#rp = 5.
#
#a = 1.
#b = 10.
#
#Q = (rm**2-rp**2+2.*rm**2*np.log(b/rm)+2.*rp**2*np.log(rp/b))/(8.*np.log(a/b))

l  = 9
rm = 0.45
rp = 0.5

a = 0.
b = 1.
Q = 0.

omq = -4.*Q/rm**2
omd = 1./2.
omr = omq / omd

bl = l*(1.-(rm/rp)**2+omr*(1.+(rm/rp)**2))*(1.-(a/b)**(2*l)) \
     +(1.-(rm/rp)**(2*l))*((rp/b)**(2*l)-(a/rm)**(2*l))/(1.-(a/b)**(2*l))

cl = l**2*omr*(1.-(rm/rp)**2+omr*(rm/rp)**2) \
     -l*omr*(1.-(a/rp)**(2*l))*(1.-(rp/b)**(2*l))/(1.-(a/b)**(2*l)) \
     +l*(1.-(rm/rp)**2+omr*(rm/rp)**2)*(1.-(rm/b)**(2*l))*(1.-(a/rm)**(2*l))/(1.-(a/b)**(2*l)) \
     -(1.-(rp/b)**(2*l))*(1.-(a/rm)**(2*l))*(1.-(rm/rp)**(2*l))/(1.-(a/b)**(2*l))

coeff = [ 1./omd**2, -bl/omd, cl ]
roots = np.roots( coeff )

print( "\n Azimuthal mode number: l =", l )
print( "\n rm/b =", rm/b )
print(   " rp/b =", rp/b )
print( "\n Growth rates:" )
print( "\n Im(omega_1) =", roots[0].imag )
print(   " Im(omega_2) =", roots[1].imag )
print( "\n Re(omega_1) =", roots[0].real )
print(   " Re(omega_2) =", roots[1].real )

## smoothing
#
#import matplotlib.pyplot as plt
#
#r = np.linspace( 0., 1., 10000 )
#
#def y_func( r ):
#    if ( rm <= r and r <= rp ):
#       return 1.
#    else:
#       return 0.
#
#y_func_vect = np.vectorize( y_func )
#
#y = y_func_vect( r )
#
#h = (rp+rm)/2.
#s = np.exp( -( (r-h)/0.025 )**50 )
#
#fg = plt.figure()
#ax = fg.add_subplot(111)
#ax.plot( r, y  )
#ax.plot( r, s )
#ax.grid()
#fg.tight_layout()
#fg.show()
