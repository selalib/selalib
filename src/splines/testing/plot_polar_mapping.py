import h5py
import numpy as np
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
# Target mapping
#-------------------------------------------------------------------------------

f0 = h5py.File('mapping_analytical_target_mesh.h5',mode='r')
f1 = h5py.File('mapping_iga_target_mesh.h5',mode='r')
f2 = h5py.File('mapping_iga_target_data.h5',mode='r')

# Mesh (from analytical mapping)
x1_target = f0['x1'].value
x2_target = f0['x2'].value

# Mesh (from iga mapping)
x1_target_iga = f1['x1'].value
x2_target_iga = f1['x2'].value

# Control points
c_x1_target = f2['c_x1'].value
c_x2_target = f2['c_x2'].value

f0.close()
f1.close()
f2.close()

#-------------------------------------------------------------------------------
# Czarny mapping
#-------------------------------------------------------------------------------

f0 = h5py.File('mapping_analytical_czarny_mesh.h5',mode='r')
f1 = h5py.File('mapping_iga_czarny_mesh.h5',mode='r')
f2 = h5py.File('mapping_iga_czarny_data.h5',mode='r')

# Mesh (from analytical mapping)
x1_czarny = f0['x1'].value
x2_czarny = f0['x2'].value

# Mesh (from iga mapping)
x1_czarny_iga = f1['x1'].value
x2_czarny_iga = f1['x2'].value

# Control points
c_x1_czarny = f2['c_x1'].value
c_x2_czarny = f2['c_x2'].value

f0.close()
f1.close()
f2.close()

#-------------------------------------------------------------------------------
# PLOTS
#-------------------------------------------------------------------------------

# Target mapping

fig = plt.figure()
ax  = fig.add_subplot(111)

# plot analytical mesh
ax.plot( x1_target, x2_target, color='b' )
ax.plot( x1_target.transpose(), x2_target.transpose(), color='b' )
# plot discrete mesh
ax.plot( x1_target_iga, x2_target_iga, '.', color='r' )
ax.plot( x1_target_iga.transpose(), x2_target_iga.transpose(), '.', color='r' )
# plot control points
ax.plot( c_x1_target.ravel(), c_x2_target.ravel(), 'x', color='k' )

plt.title('Target mapping: mesh and control points')

fig.show()

# Czarny mapping

fig = plt.figure()
ax  = fig.add_subplot(111)

# plot analytical mesh
ax.plot( x1_czarny, x2_czarny, color='b' )
ax.plot( x1_czarny.transpose(), x2_czarny.transpose(), color='b' )
# plot discrete mesh
ax.plot( x1_czarny_iga, x2_czarny_iga, '.', color='r' )
ax.plot( x1_czarny_iga.transpose(), x2_czarny_iga.transpose(), '.', color='r' )
# plot control points
ax.plot( c_x1_czarny.ravel(), c_x2_czarny.ravel(), 'x', color='k' )

plt.title('Czarny mapping: mesh and control points')

fig.show()
