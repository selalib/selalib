Installation
************

Selalib is available as `an archive file </releases/selalib-0.5.0.tar.gz>`_

Install Dependencies ::

 $ apt-get install gfortran openmpi-bin libhdf5-openmpi-dev \
                  libopenmpi-dev doxygen liblapack-dev \
                  libfftw3-dev cmake gnuplot

Build the library and run examples::
       
 $ cmake .
 $ make 
 $ make examples
 $ ./examples/landau | gnuplot

