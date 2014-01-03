#guiding center simulations

ctest --verbose -R sim2d_gc_cart
#./bin/test_2d_gc_cartesian ../selalib/prototype/src/simulation/gcsim2d_cartesian_input
gnuplot -persist ../selalib/prototype/src/simulation/gcsim2d_cartesian.gnu

ctest --verbose -R sim2d_gc_polar
#./bin/test_2d_gc_polar ../selalib/prototype/src/simulation/gcsim2d_polar_input
gnuplot -persist ../selalib/prototype/src/simulation/gcsim2d_polar.gnu

./bin/test_2d_gc_curvilinear ../selalib/prototype/src/simulation/gcsim2d_curvilinear_input
gnuplot -persist ../selalib/prototype/src/simulation/gcsim2d_curvilinear.gnu


#vlasov poisson 2D for keen and high order splitting (MPI)
ctest --verbose -R sim2d_vp_cart
#./bin/test_2d_vp_cartesian ../selalib/prototype/src/simulation/vpsim2d_cartesian_input
gnuplot -persist ../selalib/prototype/src/simulation/vpsim2d_cartesian.gnu

#for keen
./bin/test_2d_vp_cartesian ../selalib/prototype/src/simulation/vpsim2d_cartesian_keen
#mpirun -np 8 ./bin/test_2d_vp_cartesian ../selalib/prototype/src/simulation/vpsim2d_cartesian_keen
gnuplot -persist ../selalib/prototype/src/simulation/vpsim2d_cartesian_keen.gnu

#for beam
./bin/test_2d_vp_cartesian ../selalib/prototype/src/simulation/vpsim2d_cartesian_beam
#mpirun -np 8 ./bin/test_2d_vp_cartesian ../selalib/prototype/src/simulation/vpsim2d_cartesian_beam
gnuplot -persist ../selalib/prototype/src/simulation/vpsim2d_cartesian_beam.gnu

#for vlasov without splitting
./bin/test_2d_vp_no_split ../selalib/prototype/src/simulation/vpsim2d_no_split_beam
gnuplot -persist ../selalib/prototype/src/simulation/vpsim2d_no_split_beam.gnu

#vlasov poisson 4D parallel for vlasov, sequential for poisson
ctest --verbose -R sim4d_vp_cart
#time mpirun -np 8 ./bin/test_4d_vp_cartesian ../selalib/prototype/src/simulation/vpsim4d_cartesian_input
gnuplot -persist ../selalib/prototype/src/simulation/vpsim4d_cartesian.gnu


#drift kinetic 4D polar
ctest --verbose -R sim4d_DK_polar
#time mpirun -np 16 ./bin/test_4d_dk_polar ../selalib/prototype/src/simulation/dksim4d_polar_input.nml
gnuplot -persist ../selalib/prototype/src/simulation/dksim4d_polar.gnu 
