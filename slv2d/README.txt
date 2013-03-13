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


