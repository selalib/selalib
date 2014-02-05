SLL_ROOT = /Users/navaro/Codes/selalib/usr
F90 = gfortran
OPT = -O3 #-xHost -ip -openmp

HDF5_INCDIR = /usr/local/include
HDF5_LIBDIR = /usr/local/lib
HDF5_LIBRARIES = -L${HDF5_LIBDIR} -lhdf5_fortran -lhdf5

MPIF90 = mpif90
   
F90FLAGS = -I${SLL_ROOT}/include/selalib -ffree-line-length-none
LDFLAGS = ${F90FLAGS} -L${SLL_ROOT}/lib 
LIBS = -lselalib -ldfftpack ${HDF5_LIBRARIES}

PROGS = landau polar_advection
all: ${PROGS}

landau: landau.F90
	$(F90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

polar_advection: polar_advection.F90
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

run:
	./landau | gnuplot
 
clean:
	rm -f $(PROG) *.o *.mod *.dat *.gnu

.SUFFIXES: $(SUFFIXES) .F90

.F90.o:
	$(F90) $(F90FLAGS) -c $<

.mod.o:
	$(F90) $(F90FLAGS) -c $*.F90
