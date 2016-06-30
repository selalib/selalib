This simulation is a Landau damping case solved with a PIC method.

- Global parameters are in module "zone".
- All functions related to particles are in the module "particules"
 - Initialization
 - Interpolation
 - Pusher
 - Deposition

to run the code in build directory
- Copy the data file

$ selalib/simulations/parallel/pic_vp_2d2v_cart_tau/input/data.nml

- Run with the command

$ mpirun -np 4 ./bin/sim_pic_vp_2d2v_cart_tau data.nml

ntau is fixed to 32 which is the max for processors count.
