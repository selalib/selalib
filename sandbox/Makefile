SLL_ROOT = /Users/navaro/Codes/selalib/usr
F90 = gfortran
OPT = -g -Wall -ffree-line-length-none -fall-intrinsics -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow -O0

HDF5_INCDIR = /usr/local/include
HDF5_LIBDIR = /usr/local/lib
HDF5_LIBRARIES = -L${HDF5_LIBDIR} -lhdf5_fortran -lhdf5

MPIF90 = mpif90
   
F90FLAGS = -DFFTW -I${SLL_ROOT}/include/selalib ${OPT}
LDFLAGS = ${F90FLAGS} 
LIBS = ${SLL_ROOT}/lib/libselalib.a \
       ${SLL_ROOT}/lib/libdfftpack.a \
       ${HDF5_LIBRARIES}

all: landau polar_advection landau_parallel parallel_advection vlasov_poisson_parallel

landau: landau.o
	$(F90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

polar_advection: polar_advection.o
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

landau_parallel: landau_parallel.o
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

parallel_advection: parallel_advection.o
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

vlasov_poisson_parallel: vlasov_poisson_parallel.o
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
