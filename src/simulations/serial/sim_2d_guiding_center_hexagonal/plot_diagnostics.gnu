# This file will create all the graphs of the diagnostics for the Guiding-center simulation
# on a hexagonal mesh. The outputs are saved on png files.

# BEFORE USING, check the following variables :
spline_deg = 2
nc_max = 160
tmax = 100

set term png;
unset key

#------------------------------------
# Time evolution files
#------------------------------------
print " TIME FILES :"
set xlabel "time t"

filename  = "./diag_gc_spline000".spline_deg."_tmax0".nc_max.".dat"


# Mass...............................
plotfile  = "Mass_evolution_time_evol.png"
set output plotfile;
titlename = "Mass evolution, for different space discretization at time tmax = ".tmax.""
print filename." -----> ".plotfile
#---
set title titlename
plot filename u 1:2 w lp;

# Minimum...............................
plotfile  = "Minimum_value_time_evol.png"
set output plotfile;
titlename = "Minimum value, for different space discretization at time tmax = ".tmax.""
print filename." -----> ".plotfile
#---
set title titlename
plot filename u 1:3 w lp

# L1 norm...............................
plotfile  = "L1norm_time_evol.png"
set output plotfile;
titlename = "L1 norm, for different space discretization at time tmax = ".tmax.""
print filename." -----> ".plotfile
#---
set title titlename
plot filename u 1:4 w lp

# L2 norm...............................
plotfile  = "L2norm_time_evol.png"
set output plotfile;
titlename = "L2 norm, for different space discretization at time tmax = ".tmax.""
print filename." -----> ".plotfile
#---
set title titlename
plot filename u 1:5 w lp

# Linf norm...............................
plotfile  = "Linfnorm_time_evol.png"
set output plotfile;
titlename = "L_inf norm, for different space discretization at time tmax = ".tmax.""
print filename." -----> ".plotfile
#---
set title titlename
plot filename u 1:6 w lp

# Energy...............................
plotfile  = "Energy_time_evol.png"
set output plotfile;
titlename = "Energy, for different space discretization at time tmax = ".tmax.""
print filename." -----> ".plotfile
#---
set title titlename
plot filename u 1:7 w lp





print ""
#------------------------------------
# Different space discretization files
#------------------------------------
print " SPACE FILES :"
set xlabel "Number of cells in an hex edge"

filename  = "./diag_gc_spline000".spline_deg."_nc.dat"


# Mass...............................
plotfile  = "Mass_evolution_space_evol.png"
set output plotfile;
titlename = "Mass evolution, for different space discretization at time tmax = ".tmax.""
print filename." -----> ".plotfile
#---
set title titlename
plot filename u 1:3 w lp;

# Minimum...............................
plotfile  = "Minimum_value_space_evol.png"
set output plotfile;
titlename = "Minimum value, for different space discretization at time tmax = ".tmax.""
print filename." -----> ".plotfile
#---
set title titlename
plot filename u 1:4 w lp

# L1 norm...............................
plotfile  = "L1norm_space_evol.png"
set output plotfile;
titlename = "L1 norm, for different space discretization at time tmax = ".tmax.""
print filename." -----> ".plotfile
#---
set title titlename
plot filename u 1:5 w lp

# L2 norm...............................
plotfile  = "L2norm_space_evol.png"
set output plotfile;
titlename = "L2 norm, for different space discretization at time tmax = ".tmax.""
print filename." -----> ".plotfile
#---
set title titlename
plot filename u 1:6 w lp

# Linf norm...............................
plotfile  = "Linfnorm_space_evol.png"
set output plotfile;
titlename = "L_inf norm, for different space discretization at time tmax = ".tmax.""
print filename." -----> ".plotfile
#---
set title titlename
plot filename u 1:7 w lp

# Energy...............................
plotfile  = "Energy_space_evol.png"
set output plotfile;
titlename = "Energy, for different space discretization at time tmax = ".tmax.""
print filename." -----> ".plotfile
#---
set title titlename
plot filename u 1:8 w lp

