SLL_ROOT = /Users/navaro/Codes/selalib/usr
F90 = gfortran
OPT = -g -Wall -ffree-line-length-none -fall-intrinsics -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow -O0

HDF5_INCDIR = /usr/local/include
HDF5_LIBDIR = /usr/local/lib
HDF5_LIBRARIES = -L${HDF5_LIBDIR} -lhdf5_fortran -lhdf5

MPIF90 = mpif90
   
F90FLAGS = -I${SLL_ROOT}/include/selalib ${OPT}
LDFLAGS = ${F90FLAGS} 
LIBS = ${SLL_ROOT}/lib/libselalib.a \
       ${SLL_ROOT}/lib/libdfftpack.a \
       ${HDF5_LIBRARIES}

PROGS = landau polar_advection landau_parallel parallel_advection
all: ${PROGS}

landau: landau.F90
	$(F90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

polar_advection: polar_advection.F90
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

landau_parallel: landau_parallel.F90
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

parallel_advection: parallel_advection.F90
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

run:
	./landau | gnuplot
 
clean:
	rm -f *.exe *.o *.mod *.dat *.gnu

.SUFFIXES: $(SUFFIXES) .F90

.F90.o:
	$(F90) $(F90FLAGS) -c $<

.mod.o:
	$(F90) $(F90FLAGS) -c $*.F90
