This simulation is the code for the paper
"""
  Uniformly accurate particle-in-cell method for the long time
  two-dimensional Valsov-Poisson equation with strong magnetic field
"""
by
Nicolas Crouseilles, Mohammed Lemou, Florian Mehats and Xiaofei Zhao

- Global parameters are in module "zone".
- All functions related to particles are in the module "particules"
	- Initialization
	- Interpolation
	- Pusher
	- Deposition

to run the code in build directory
- Copy the data file 

$ selalib/simulations/parallel/pic_vp_2d2v_cart_tau/data.nml

- Run with the command

$ mpirun -np 4 ./bin/sim_pic_vp_2d2v_cart_tau data.nml

ntau value is the max for processors count.

You get the first fourier modes of energy in file energy.dat
