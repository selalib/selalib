import h5py
import numpy as np
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
# Circle mapping
#-------------------------------------------------------------------------------

f0 = h5py.File('mapping_analytical_circle.h5',mode='r')
f1 = h5py.File('mapping_iga_circle.h5',mode='r')

# Load data from analytical mapping
x1_circle_analytical       = f0['x1'].value
x2_circle_analytical       = f0['x2'].value
jacobian_circle_analytical = f0['jacobian'].value

# Load data from discrete IGA mapping
x1_circle_iga       = f1['x1'].value
x2_circle_iga       = f1['x2'].value
c_x1_circle_iga     = f1['c_x1'].value
c_x2_circle_iga     = f1['c_x2'].value
jacobian_circle_iga = f1['jacobian'].value

f0.close()
f1.close()

#-------------------------------------------------------------------------------
# Target mapping
#-------------------------------------------------------------------------------

f0 = h5py.File('mapping_analytical_target.h5',mode='r')
f1 = h5py.File('mapping_iga_target.h5',mode='r')

# Load data from analytical mapping
x1_target_analytical       = f0['x1'].value
x2_target_analytical       = f0['x2'].value
jacobian_target_analytical = f0['jacobian'].value

# Load data from discrete IGA mapping
x1_target_iga       = f1['x1'].value
x2_target_iga       = f1['x2'].value
c_x1_target_iga     = f1['c_x1'].value
c_x2_target_iga     = f1['c_x2'].value
jacobian_target_iga = f1['jacobian'].value

f0.close()
f1.close()

#-------------------------------------------------------------------------------
# Czarny mapping
#-------------------------------------------------------------------------------

f0 = h5py.File('mapping_analytical_czarny.h5',mode='r')
f1 = h5py.File('mapping_iga_czarny.h5',mode='r')

# Load data from analytical mapping
x1_czarny_analytical       = f0['x1'].value
x2_czarny_analytical       = f0['x2'].value
jacobian_czarny_analytical = f0['jacobian'].value

# Load data from discrete IGA mapping
x1_czarny_iga       = f1['x1'].value
x2_czarny_iga       = f1['x2'].value
c_x1_czarny_iga     = f1['c_x1'].value
c_x2_czarny_iga     = f1['c_x2'].value
jacobian_czarny_iga = f1['jacobian'].value

f0.close()
f1.close()

#-------------------------------------------------------------------------------
# OUTPUT
#-------------------------------------------------------------------------------

print()
print( ' Circle mapping' )
print( ' ==============' )
print()
print(' Maximum absolute errors between analytical and discrete mapping:')
print()
print( ' x:', np.amax( np.abs( x1_circle_analytical - x1_circle_iga ) ) )
print( ' y:', np.amax( np.abs( x2_circle_analytical - x2_circle_iga ) ) )
print( ' J:', np.amax( np.abs( jacobian_circle_analytical - jacobian_circle_iga ) ) )

print()
print( ' Target mapping' )
print( ' ==============' )
print()
print(' Maximum absolute errors between analytical and discrete mapping:')
print()
print( ' x:', np.amax( np.abs( x1_target_analytical - x1_target_iga ) ) )
print( ' y:', np.amax( np.abs( x2_target_analytical - x2_target_iga ) ) )
print( ' J:', np.amax( np.abs( jacobian_target_analytical - jacobian_target_iga ) ) )

print()
print( ' Czarny mapping' )
print( ' ==============' )
print()
print(' Maximum absolute errors between analytical and discrete mapping:')
print()
print( ' x:', np.amax( np.abs( x1_czarny_analytical - x1_czarny_iga ) ) )
print( ' y:', np.amax( np.abs( x2_czarny_analytical - x2_czarny_iga ) ) )
print( ' J:', np.amax( np.abs( jacobian_czarny_analytical - jacobian_czarny_iga ) ) )

#-------------------------------------------------------------------------------
# PLOTS
#-------------------------------------------------------------------------------

# Circle mapping

fig = plt.figure()
ax  = fig.add_subplot(111)

# plot analytical mesh
ax.plot( x1_circle_analytical, x2_circle_analytical, color='b' )
ax.plot( x1_circle_analytical.transpose(), x2_circle_analytical.transpose(), color='b' )
# plot discrete mesh
ax.plot( x1_circle_iga, x2_circle_iga, '.', color='r' )
ax.plot( x1_circle_iga.transpose(), x2_circle_iga.transpose(), '.', color='r' )
# plot control points
ax.plot( c_x1_circle_iga.ravel(), c_x2_circle_iga.ravel(), 'x', color='k' )

plt.title('Circle mapping: mesh and control points')

fig.show()

# Target mapping

fig = plt.figure()
ax  = fig.add_subplot(111)

# plot analytical mesh
ax.plot( x1_target_analytical, x2_target_analytical, color='b' )
ax.plot( x1_target_analytical.transpose(), x2_target_analytical.transpose(), color='b' )
# plot discrete mesh
ax.plot( x1_target_iga, x2_target_iga, '.', color='r' )
ax.plot( x1_target_iga.transpose(), x2_target_iga.transpose(), '.', color='r' )
# plot control points
ax.plot( c_x1_target_iga.ravel(), c_x2_target_iga.ravel(), 'x', color='k' )

plt.title('Target mapping: mesh and control points')

fig.show()

# Czarny mapping

fig = plt.figure()
ax  = fig.add_subplot(111)

# plot analytical mesh
ax.plot( x1_czarny_analytical, x2_czarny_analytical, color='b' )
ax.plot( x1_czarny_analytical.transpose(), x2_czarny_analytical.transpose(), color='b' )
# plot discrete mesh
ax.plot( x1_czarny_iga, x2_czarny_iga, '.', color='r' )
ax.plot( x1_czarny_iga.transpose(), x2_czarny_iga.transpose(), '.', color='r' )
# plot control points
ax.plot( c_x1_czarny_iga.ravel(), c_x2_czarny_iga.ravel(), 'x', color='k' )

plt.title('Czarny mapping: mesh and control points')

fig.show()
