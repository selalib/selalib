SRCS_DIR = /Users/navaro/Codes/selalib/prototype/src/package/examples
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

all: landau                  \
     polar_advection         \
     landau_parallel         \
     parallel_advection      \
     poisson_2d_parallel     \
     vlasov_poisson_parallel \
     vlasov_maxwell_parallel 

landau: ${SRCS_DIR}/landau.F90
	$(F90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

polar_advection: ${SRCS_DIR}/polar_advection.F90
	$(F90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

landau_parallel: ${SRCS_DIR}/landau_parallel.F90
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

parallel_advection: ${SRCS_DIR}/parallel_advection.F90
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

vlasov_poisson_parallel: ${SRCS_DIR}/vlasov_poisson_parallel.F90
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)
 
vlasov_maxwell_parallel: ${SRCS_DIR}/vlasov_maxwell_parallel.F90
	$(MPIF90) -DFFTW $(LDFLAGS) -o $@.exe $^ $(LIBS)

poisson_2d_parallel: ${SRCS_DIR}/poisson_2d_parallel.F90
	$(MPIF90) $(LDFLAGS) -o $@.exe $^ $(LIBS)

clean:
	rm -rf *.exe *.o *.mod *.dat *.gnu *.exe.dSYM 00* fort.*
