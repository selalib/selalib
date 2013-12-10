import numpy as np
import matplotlib.pylab as plb
import matplotlib.pyplot as plt
import utils_func as ut 

CGresu = plb.mlab.load('thdiag.dat',skiprows=4)

#--> time evolution of the mode 3
CGtime     = CGresu[:,0]
real_mode3 = CGresu[:,12]
imag_mode3 = CGresu[:,13]
abs2_mode3 = real_mode3**2 + imag_mode3**2

#--> analytical frequency and growth rate of the mode 3
real_mode3_anal = -0.5456
imag_mode3_anal = 0.2267
val0            = abs2_mode3[1]
abs2_mode3_anal = val0*np.exp(0.2267*CGtime)

#--> comparison simulation/analytical growth rate
plt.figure()
plt.hold(True)
plt.grid('on')
str_leg = []
plt.semilogy(CGtime,abs2_mode3,'b')
str_leg.append('mode3 num')
plt.semilogy(CGtime,abs2_mode3_anal,'g')
str_leg.append('anal growth rate')

#--> numerical computation of the growth rate
[t_init,t_end] = ut.Ask_time(CGtime,ask_diagstep=False,msg='for growth rate computation')
[p,s] = np.polyfit(CGtime[t_init:t_end+1], \
         np.log(abs2_mode3[t_init:t_end+1]),1)
growthrate_mode3 = p
val0_mode3       = s
print 'growthrate_mode3 = ' + str(p)
print 'val0_mode3 = ' + str(s)
plt.semilogy(CGtime,np.exp(val0_mode3) * \
             np.exp(growthrate_mode3*CGtime),'r--')
str_leg.append('num growth rate')
plt.legend(str_leg)
plt.title('num growth rate = ' + str(growthrate_mode3) + \
          ' compared to ' + str(imag_mode3_anal))
plt.hold(False)
