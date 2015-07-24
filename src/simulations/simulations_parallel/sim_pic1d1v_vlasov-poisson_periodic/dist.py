from pylab import *
import numpy as np
import matplotlib.pyplot as plt
y, x = np.meshgrid(np.linspace(0, 13,300), np.linspace(-4.5, 4.5,200))
ax = plt.subplot(111)
xlabel('position')
ylabel('velocity')
j=1
string='distribution__'
file_number=10000+j
file_number=str(file_number)
file_number=file_number[1:5]
file_string=string+file_number+'.xmf'
f=np.genfromtxt(file_string,skip_header=16,skip_footer=5)
m=[f[0:300]]
for i in range(2,201):
     m.append(f[(i-1)*300:300*i])

mm1=np.asarray(m)

quad = plt.pcolormesh(y,x,mm1)
plt.colorbar()
plt.ion()


for j in range(1,200):
  file_number=10000+j
  file_number=str(file_number)
  file_number=file_number[1:5]
  file_string=string+file_number+'.xmf'
  f=np.genfromtxt(file_string,skip_header=16,skip_footer=5)
  m=[f[0:300]]
  for i in range(2,201):
      m.append(f[(i-1)*300:300*i])

  mm=np.asarray(m)
  quad.set_array(mm.ravel())
 # pcolormesh(y,x,mm)

  plt.title('bum on tail'+file_number)
  plt.draw()
 
plt.show()
plt.ioff()



