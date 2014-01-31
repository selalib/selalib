SLL_ROOT = /Users/navaro/Codes/selalib/usr
F90 = gfortran
OPT = -g -Wall -ffree-line-length-none -fall-intrinsics -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow -O0

OPT = -O3 -ffree-line-length-none 

HDF5_INCDIR = /usr/local/include
HDF5_LIBDIR = /usr/local/lib
HDF5_LIBRARIES = -L${HDF5_LIBDIR} -lhdf5_fortran -lhdf5

FFTW_INCDIR = /usr/local/include
FFTW_LIBDIR = /usr/local/lib
FFTW_LIBRARIES = -L${FFTW_LIBDIR} -lfftw3 

MPIF90 = mpif90
   
F90FLAGS = -I${SLL_ROOT}/include/selalib ${OPT}
LDFLAGS = ${F90FLAGS} 
LIBS = ${SLL_ROOT}/lib/libselalib.a \
       ${SLL_ROOT}/lib/libdfftpack.a \
       ${HDF5_LIBRARIES} ${FFTW_LIBRARIES}

all: landau polar_advection landau_parallel \
     parallel_advection vlasov_poisson_parallel \
     poisson_2d_parallel

landau: landau.o
	$(F90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

polar_advection: polar_advection.o
	$(F90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

landau_parallel: landau_parallel.o
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

parallel_advection: parallel_advection.o
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

vlasov_poisson_parallel: vlasov_poisson_parallel.o
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)
 
vlasov_maxwell_parallel: vlasov_maxwell_parallel.o
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

poisson_2d_parallel: poisson_2d_parallel.o
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

clean:
	rm -f *.exe *.o *.mod *.dat *.gnu

.SUFFIXES: $(SUFFIXES) .F90

.F90.o:
	$(F90) $(F90FLAGS) -c $<

.mod.o:
	$(F90) $(F90FLAGS) -c $*.F90

landau.o: landau.F90
landau_parallel.o: landau_parallel.F90
parallel_advection.o: parallel_advection.F90
polar_advection.o: polar_advection.F90
vlasov_poisson_parallel.o: vlasov_poisson_parallel.F90
vlasov_maxwell_parallel.o: vlasov_maxwell_parallel.F90
poisson_2d_parallel.o : poisson_2d_parallel.F90
