# coding: utf-8
from jorek.plot.quadrangles import quad_plot
from matplotlib import pyplot as plt
from matplotlib.tri import UniformTriRefiner
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import sys

# -----------------------------------------------
try:
    filename = sys.argv[1]
except:
    filename = "bin/poisson_2d_0.tp"
# -----------------------------------------------

quad = quad_plot(filename=filename)

x = quad.triang.x
y = quad.triang.y

uh =quad.variables[0]
eh =quad.variables[1]

plt.tripcolor(quad.triang, uh, shading='gouraud', cmap=plt.cm.rainbow); plt.colorbar()
plt.title("$u_h$")
plt.savefig("uh.png")
plt.clf()

plt.tripcolor(quad.triang, eh, shading='gouraud', cmap=plt.cm.rainbow); plt.colorbar()
plt.title("$e_h=u-u_h$")
plt.savefig("eh.png")
plt.clf()

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(x, y, uh, cmap=cm.jet, linewidth=0.2)
plt.title("$u_h$")
plt.savefig("uh-3d.png")
plt.clf()

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(x, y, eh, cmap=cm.jet, linewidth=0.2)
plt.title("$e_h$")
plt.savefig("eh-3d.png")
plt.clf()


#plt.clf()

#triang = quad.triang
#values = quad.variables[0]
#
#refiner = UniformTriRefiner(triang)
#triang_new, values= refiner.refine_field(values, subdiv=3)
#
#values = np.array(values)
#print values.min(), values.max()
#
#plt.tripcolor(triang_new, values, shading='gouraud', cmap=plt.cm.rainbow)
#plt.savefig("uh_new.png")
#
