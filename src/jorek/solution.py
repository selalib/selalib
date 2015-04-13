from sympy import *

label = "lr_" # local variable

R   = Symbol(label+"R")
Z   = Symbol(label+"Z")

a      = Symbol(label+"a")
acenter   = Symbol(label+"acenter")
R0  = Symbol(label+"R0")

K  = eye(2)

#u = 1 - ((R-R0)**2 + Z**2)/a**2
u = (1 - ((R-R0)**2 + Z**2)/a**2)*(((R-R0)**2 + Z**2)- acenter**2)

u_R   = diff(u,R)
u_Z   = diff(u,Z)

grad_u = Matrix([u_R, u_Z])

KU =K.multiply(grad_u)

div_KU = diff(KU[0],R)  + diff(KU[1],Z)

f = - div_KU

print "u"
print u
print "====="
print "grad u"
print grad_u
print "====="
print "K"
print K
print "====="
print "KU"
print KU
print "====="
print "====="
print "====="
print "f"
print f
