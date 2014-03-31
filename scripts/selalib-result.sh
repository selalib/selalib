# swirling deformation flow (analytical field)
./bin/test_2d_af_cartesian ../selalib/prototype/src/simulation/afsim2d_cartesian_sdf
gnuplot -persist ../selalib/prototype/src/simulation/afsim2d_cartesian_sdf.gnu

#../selalib/prototype/src/simulation/simulation_2d_analytic_field_curvilinear.F90
#./bin/test_2d_af_curvilinear ../selalib/prototype/src/simulation/afsim2d_curvilinear_sdf


#guiding center simulations

ctest --verbose -R sim2d_gc_cart
#./bin/test_2d_gc_cartesian ../selalib/prototype/src/simulation/gcsim2d_cartesian_input
gnuplot -persist ../selalib/prototype/src/simulation/gcsim2d_cartesian.gnu

ctest --verbose -R sim2d_gc_polar
#./bin/test_2d_gc_polar ../selalib/prototype/src/simulation/gcsim2d_polar_input
gnuplot -persist ../selalib/prototype/src/simulation/gcsim2d_polar.gnu

#curvilinear colella
#./bin/test_2d_gc_curvilinear ../selalib/prototype/src/simulation/gcsim2d_curvilinear_input
#gnuplot -persist ../selalib/prototype/src/simulation/gcsim2d_curvilinear.gnu

#curvilinear cartesian
#./bin/test_2d_gc_curvilinear ../selalib/prototype/src/simulation/gcsim2d_curvilinear_input_cart
#gnuplot -persist ../selalib/prototype/src/simulation/gcsim2d_curvilinear_cart.gnu

#curvilinear polar
#./bin/test_2d_gc_curvilinear ../selalib/prototype/src/simulation/gcsim2d_curvilinear_input_polar
#gnuplot -persist ../selalib/prototype/src/simulation/gcsim2d_curvilinear_polar.gnu



#vlasov poisson 2D for keen and high order splitting (MPI)
ctest --verbose -R sim2d_vp_cart
#./bin/test_2d_vp_cartesian ../selalib/prototype/src/simulation/vpsim2d_cartesian_input
gnuplot -persist ../selalib/prototype/src/simulation/vpsim2d_cartesian.gnu

#for bump on tail
mpirun -np 2 ./bin/test_2d_vp_cartesian ../selalib/prototype/src/simulation/vpsim2d_cartesian_bot
gnuplot -persist ../selalib/prototype/src/simulation/vpsim2d_cartesian_bot.gnu

#for two stream instability
mpirun -np 2 ./bin/test_2d_vp_cartesian ../selalib/prototype/src/simulation/vpsim2d_cartesian_tsi
gnuplot -persist ../selalib/prototype/src/simulation/vpsim2d_cartesian_tsi.gnu


#for keen
./bin/test_2d_vp_cartesian ../selalib/prototype/src/simulation/vpsim2d_cartesian_keen
#mpirun -np 8 ./bin/test_2d_vp_cartesian ../selalib/prototype/src/simulation/vpsim2d_cartesian_keen
gnuplot -persist ../selalib/prototype/src/simulation/vpsim2d_cartesian_keen.gnu

#for beam
./bin/test_2d_vp_cartesian ../selalib/prototype/src/simulation/vpsim2d_cartesian_beam
#mpirun -np 8 ./bin/test_2d_vp_cartesian ../selalib/prototype/src/simulation/vpsim2d_cartesian_beam
gnuplot -persist ../selalib/prototype/src/simulation/vpsim2d_cartesian_beam.gnu

#for vlasov without splitting
ctest --verbose -R sim2d_vp_no_split
#./bin/test_2d_vp_no_split ../selalib/prototype/src/simulation/vpsim2d_no_split_beam
gnuplot -persist ../selalib/prototype/src/simulation/vpsim2d_no_split_beam.gnu

#vlasov poisson 4D parallel for vlasov, sequential for poisson
ctest --verbose -R sim4d_vp_cart
#time mpirun -np 8 ./bin/test_4d_vp_cartesian ../selalib/prototype/src/simulation/vpsim4d_cartesian_input
gnuplot -persist ../selalib/prototype/src/simulation/vpsim4d_cartesian.gnu


#drift kinetic 4D polar
ctest --verbose -R sim4d_DK_polar
#time mpirun -np 16 ./bin/test_4d_dk_polar ../selalib/prototype/src/simulation/dksim4d_polar_input.nml
gnuplot -persist ../selalib/prototype/src/simulation/dksim4d_polar.gnu 

#drift kinetic 4D polar one mu
ctest --verbose -R sim4d_DK_polar_one_mu
#time mpirun -np 16 ./bin/test_4d_dk_polar_one_mu ../selalib/prototype/src/simulation/dksim4d_polar_one_mu.nml
gnuplot -persist ../selalib/prototype/src/simulation/dksim4d_polar_one_mu.gnu 

#drift kinetic 4D polar one mu
ctest --verbose -R sim4d_DK_polar_multi_mu
#time mpirun -np 16 ./bin/test_4d_dk_polar_multi_mu ../selalib/prototype/src/simulation/dksim4d_polar_multi_mu.nml
gnuplot -persist ../selalib/prototype/src/simulation/dksim4d_polar_multi_mu.gnu 




