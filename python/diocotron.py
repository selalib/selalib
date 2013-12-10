from sympy.solvers import solve
from sympy import Symbol
import numpy as np

ell    = 3
rm     = 4.
rp     = 6.
omegad = 0.5
a      = 1
b      = 10

bell = ell*(1-(rm/rp)**2)*(1-(a/b)**(2*ell))
bell = bell+(1-(rm/rp)**(2*ell)) * \
    ((rp/b)**(2*ell)-(a/rm)**(2*ell))*(1-(a/b)**(2*ell))**(-1)
cell = ell*(1-(rm/rp)**2)*(1-(rm/b)**(2*ell))*(1-(a/rm)**(2*ell))
cell = cell-(1-(rp/b)**(2*ell))*(1-(a/rm)**(2*ell)) * \
    (1-(rm/rp)**(2*ell))*(1-(a/b)**(2*ell))**(-1)

x  = Symbol('x')
A  = solve((x/omegad)**2-bell*(x/omegad)+cell,x)
gr = 2*np.array(A)

print 'gr = ' + str(gr)
