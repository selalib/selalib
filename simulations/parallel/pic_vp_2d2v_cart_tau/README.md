UNIFORMLY ACCURATE PARTICLE-IN-CELL METHOD FOR HIGHLY OSCILLATORY VLASOV-POISSON EQUATION
=========================================================================================

This simulation is the code for the paper
"""
  Uniformly accurate particle-in-cell method for the long time
  two-dimensional Valsov-Poisson equation with strong magnetic field
"""
by
Nicolas Crouseilles, Mohammed Lemou, Florian Mehats and Xiaofei Zhao

This application is running the numerical simulation of a Vlasov-Poisson
equation modeling charged particles in a beam submitted to a highly
oscillatory external electric field. The numerical scheme is designed
to be uniformly accurate with respect to
the size of the fast time oscillations of the solution, which means
that no time step refinement is required to simulate the problem.

The scheme combines the Particle-In-Cell method with a class
of Uniformly Accurate (UA) time integrators to solve the characteristics.
These UA time integrators are derived by means of a two-scale
formulation of the characteristics, with the introduction of an
additional periodic variable.


to run the code in build directory
- Copy the data file 

$ selalib/simulations/parallel/pic_vp_2d2v_cart_tau/data.nml

or create one following instructions below.

- Run with the command

$ mpirun -np 4 ./bin/sim_pic_vp_2d2v_cart_tau data.nml

ntau value is the max for processors count.

You get the first fourier modes of energy in file energy.dat

You can run the numerical experiments changing the following parameters:

- `nstepmax` : maximum time steps.
- `nx`       : size in x.	
- `ny`       : size in y.
- `tfinal`   : finaltime.
- `ntau`     : number of points for time discretization.
- `npm`      : number of particles by cell.
- `ep`       : epsilon.
- `alpha`    : perturbation amplitude.
- `kx`       : wave number in x.
- `ky`       : wave number in y.
- `dt`       : time step.
- `plot`     : (plot rho and df) (.true. or .false.)

The data file s a fortran namelist and can be:

$donnees
  nstepmax = 10000
  nx       = 64	
  ny       = 32
  tfinal   = 25.
  ntau     = 16
  npm      = 100
  ep       = 0.005
  alpha    = 0.02
  kx       = 0.5
  ky       = 1.0
  dt       = 0.1
  plot     = .false.
$end
