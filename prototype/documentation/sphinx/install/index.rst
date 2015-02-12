Build the library and run examples
**********************************

Selalib is available as `an archive file </releases/selalib-0.5.0.tar.gz>`_

Install Dependencies ::

 $ apt-get install gfortran openmpi-bin libhdf5-openmpi-dev \
                  libopenmpi-dev doxygen liblapack-dev \
                  libfftw3-dev cmake gnuplot

Build the library and run examples::
       
 $ cmake .
 $ make 
 $ make examples

Landau damping::

 $ ./examples/landau | gnuplot
 $ mpirun -np 4 ./examples/landau_parallel | gnuplot

Simple advection::

 $ ./examples/polar_advection
 $ gnuplot f_polar_advection.gnu 
 $ mpirun -np 4 ./examples/parallel_advection
 $ gnuplot f_parallel.gnu 

Vlasov-Poisson in 4d space ::

 $ mpirun -np 4 ./examples/vlasov_poisson_parallel | gnuplot
