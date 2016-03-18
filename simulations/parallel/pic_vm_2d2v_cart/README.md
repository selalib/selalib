This simulation is a Landau damping case solved with a PIC method.

- Global parameters are in module "zone".
- All functions related to particles are in the module "particules"
 - Initialization
 - Interpolation
 - Pusher
 - Deposition

to run the code in build directory
- Copy the data file

$ selalib/simulations/parallel/pic_vm_2d2v_cart/input/data.nml

- Run with the command

$ ./bin/sim_pic_vm_2d2v_cart data.nml

