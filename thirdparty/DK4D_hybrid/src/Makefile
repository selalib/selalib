EXEC = DK4D_hybrid.exe

SRCS_DIR = ${PWD}
SLL_ROOT = ${HOME}/SELALIB_proto_lib

#OPT = -g -ffree-line-length-none
OPT = -cpp -g -check bounds -fpe0

#VG#HDF5_INCDIR = -I/opt/local/include 
#VG#HDF5_LIBDIR = /opt/local/lib 
HDF5_INCDIR = -I/gpfslocal/pub/hdf5/build_1.8.10_intel_intelmpi/include
HDF5_LIBDIR = /gpfslocal/pub/hdf5/build_1.8.10_intel_intelmpi/lib
HDF5_LIBRARIES = -L${HDF5_LIBDIR} -lhdf5_fortran -lhdf5

ZLIB_LIBRARY = -lz
#VG#ZLIB_LIBRARY = /usr/lib/libz.so

#FFTW_INCDIR = /usr/local/include
#FFTW_LIBDIR = /usr/local/lib
#FFTW_LIBRARIES = -L${FFTW_LIBDIR} -lfftw3 

MPIF90 = mpif90

OBJS = utils_DK4D.o input_DK4D.o mesh_DK4D.o magnetconf_DK4D.o equilibrium_DK4D.o fdistribu_DK4D.o field3d_DK4D.o elecpot_DK4D_hybrid.o elecpot_DK4D_polar.o QNsolver_DK4D_hybrid.o QNsolver_DK4D_polar.o vlasov_DK4D_hybrid.o vlasov_DK4D_polar.o diagnostics_DK4D.o simulation_DK4D.o main_DK4D.o
   
F90FLAGS = ${OPT} -I${SLL_ROOT}/include/selalib
LDFLAGS = ${F90FLAGS} 
LIBS = ${SLL_ROOT}/lib/libselalib.a \
       ${SLL_ROOT}/lib/libdfftpack.a \
       ${HDF5_LIBRARIES} ${FFTW_LIBRARIES} ${ZLIB_LIBRARY}

all: $(EXEC)

$(EXEC): $(OBJS)
	$(MPIF90) $(OBJS) $(LIBS) $(LDFLAGS) -o $@

diagnostics_DK4D.o : diagnostics_DK4D.F90
equilibrium_DK4D.o : equilibrium_DK4D.F90
elecpot_DK4D_hybrid.o : elecpot_DK4D_hybrid.F90
elecpot_DK4D_polar.o : elecpot_DK4D_polar.F90
fdistribu_DK4D.o : fdistribu_DK4D.F90
field3d_DK4D.o : field3d_DK4D.F90
input_DK4D.o : input_DK4D.F90
magnetconf_DK4D.o : magnetconf_DK4D.F90
main_DK4D.o : main_DK4D.F90
mesh_DK4D.o : mesh_DK4D.F90
QNsolver_DK4D_hybrid.o : QNsolver_DK4D_hybrid.F90
QNsolver_DK4D_polar.o : QNsolver_DK4D_polar.F90
simulation_DK4D.o : simulation_DK4D.F90
utils_DK4D.o : utils_DK4D.F90
vlasov_DK4D_hybrid.o : vlasov_DK4D_hybrid.F90
vlasov_DK4D_polar.o : vlasov_DK4D_polar.F90

$(OBJS):
	$(MPIF90) -c $(F90FLAGS) $(HDF5_INCDIR) $< -o $@

%o: %.F90
	$(MPIF90) -c $(F90FLAGS) $(HDF5_INCDIR) $< -o $@

clean:
	rm -rf *.exe *.o *.mod *.dat fort.*
