SLV2D 
=====

compilation
-----------

mkdir build
cd build
cmake <the path of this directory>/src
make 
make test

To launch a landau test case
----------------------------

mpirun -np 4 ./bin vp2d ../input/landau.nml

plot energy
-----------

gnuplot> plot 'thf.dat' u 13:12 w lines



vp4d_transpose : Vlasov-Poisson in 4D phase space. Transposition with Eric subroutines
vp4d_remapper : Vlasov-Poisson in 4D phase space. Transposition with sll_remapper
vm4d_transpose : Vlasov-Maxwell in 4D phase space. Transposition with Eric subroutines
vm4d_remapper : Vlasov-Maxwell in 4D phase space. Transposition with sll_remapper
vm4d_spectral : We use FFTs for advections.
