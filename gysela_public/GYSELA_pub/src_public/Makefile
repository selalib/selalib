#**************************************************************
#  Copyright Euratom-CEA
#  Authors : 
#     Virginie Grandgirard (virginie.grandgirard@cea.fr)
#     Chantal Passeron (chantal.passeron@cea.fr)
#     Guillaume Latu (guillaume.latu@cea.fr)
#     Xavier Garbet (xavier.garbet@cea.fr)
#     Philippe Ghendrih (philippe.ghendrih@cea.fr)
#     Yanick Sarazin (yanick.sarazin@cea.fr)
#  
#  This code GYSELA (for GYrokinetic SEmi-LAgrangian) 
#  is a 5D gyrokinetic global full-f code for simulating 
#  the plasma turbulence in a tokamak.
#  
#  This software is governed by the CeCILL-B license 
#  under French law and abiding by the rules of distribution 
#  of free software.  You can  use, modify and redistribute 
#  the software under the terms of the CeCILL-B license as 
#  circulated by CEA, CNRS and INRIA at the following URL
#  "http://www.cecill.info". 
#**************************************************************

# HDF5 librariries root
BASEHDF5=/usr/local/hdf5

# MPI fortran
F90=/opt/intel9/openmpi/bin/mpif90

# Fortran compiler options
F90FLAGS = -cpp -O3 -DCTIMER -I$(BASEHDF5)/lib

# preprocessor directive
TIMERFLAG = -DTIMER

# OpenMP flag
OMPFLAG = -openmp

# R3 FLAG
R3FLAG = -DUSE_R3_INFO_MPI

# MPICH library
LIBS = -L/opt/intel9/openmpi/lib -lmpi

# zlib libray used by HDF5 functions
LDFLAGS = ${OMPFLAG} -lz 

# C directive
CFLAGS = -O3 -DNORMAL

# HDF5 libraries
HDF5INCLUDE = -I$(BASEHDF5)/lib
HDF5LIB     = $(BASEHDF5)/lib/libhdf5_fortran.a $(BASEHDF5)/lib/libhdf5.a 
NOHDF5FLAG  = -DNOHDF5

BINDIR  = .
EXEC    = $(BINDIR)/gysela.exe
EXEC_TI = $(BINDIR)/gysela_ti.exe

OBJSOMP = $(OBJS1) $(OBJS3) $(OBJS5)

OBJS1 = prec_const.o globals.o clock.o mem_alloc.o \
	utils.o MPIutils.o OMPutils.o \
	geometry.o gys_alloc.o local_spline.o 

OBJS2 = HDF5_io.o resu3D_saving.o f5D_saving.o

OBJS3 = bessel.o LUsolving.o read_data.o spline1d.o \
	interpolation.o coord_system.o integration.o\
	read_profile.o init_profile.o init_current.o \
	init_magnetic.o \
	fft_NRF90.o fft.o gyroaverage.o \
	fequil4d.o fdistribu5d.o \
	poisson.o efield.o \
	filter.o advec1D_BSL.o \
	advec2D_BSL.o \
	diag3D.o Pcross_section.o physics.o

OBJS4 = read_write.o output_saving.o #r3_info.o

OBJS5 = diagnostics.o vlasov_solving.o gysela.o 

OBJSC = ctime.o #r3_ctlk.o

OBJS = $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJSC)
SOURCES = $(OBJS1:.o=.f90) $(OBJS2:.o=.f90) $(OBJS3:.o=.f90) $(OBJS4:.o=.f90) $(OBJS5:.o=.f90) $(BIN:=.f90)

ifeq ($(MAKECMDGOALS),noHDF5)
	OBJS = $(OBJS1) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJSC)
	SOURCES =  $(OBJS1:.o=.f90) $(OBJS3:.o=.f90) $(OBJS4:.o=.f90) $(OBJS5:.o=.f90) $(BIN:=.f90)

	NOHDF5 = $(NOHDF5FLAG)
	HDF5LIBS =
	HDF5INCLUDES =
else
	OBJS = $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJSC)
	SOURCES = $(OBJS1:.o=.f90) $(OBJS2:.o=.f90) $(OBJS3:.o=.f90) $(OBJS4:.o=.f90) $(OBJS5:.o=.f90) $(BIN:=.f90)

	NOHDF5 = 
	HDF5LIBS = $(HDF5LIB)
	HDF5INCLUDES = $(HDF5INCLUDE)

endif

.SUFFIXES:
.SUFFIXES: .F90  .f90 .o

debug:$(EXEC)


noHDF5:$(EXEC)	

ti: $(EXEC_TI)

$(EXEC): $(OBJS)
	$(F90) $(OMPFLAG) $(OBJS) $(LIBS) $(HDF5LIBS) $(LDFLAGS) -o $@ 

$(EXEC_TI): $(OBJS)
	$(F90) $(OMPFLAG) $(OBJS) $(LIBS) $(HDF5LIBS) $(LDFLAGS) -o $@ 

	rm -f $(OBJS) *.o *~ *.mod  core
          prec_const.o : prec_const.f90
             globals.o : globals.f90 prec_const.f90
               clock.o : clock.f90 globals.f90 prec_const.f90
           mem_alloc.o : mem_alloc.f90 globals.f90 prec_const.f90
               utils.o : utils.f90 mem_alloc.f90 prec_const.f90
            MPIutils.o : MPIutils.f90 globals.f90 mem_alloc.f90 prec_const.f90
            OMPutils.o : OMPutils.f90 globals.f90 mem_alloc.f90 prec_const.f90
            geometry.o : geometry.f90 globals.f90 mem_alloc.f90 prec_const.f90\
                         utils.f90
           gys_alloc.o : gys_alloc.f90 geometry.f90 globals.f90 mem_alloc.f90\
                         MPIutils.f90
        local_spline.o : local_spline.f90 geometry.f90 globals.f90 mem_alloc.f90\
                         prec_const.f90
             HDF5_io.o : HDF5_io.f90
       resu3D_saving.o : resu3D_saving.f90
          f5D_saving.o : f5D_saving.f90
              bessel.o : bessel.f90 mem_alloc.f90 prec_const.f90
           LUsolving.o : LUsolving.f90
           read_data.o : read_data.f90 globals.f90 MPIutils.f90 prec_const.f90
            spline1d.o : spline1d.f90 geometry.f90 globals.f90 mem_alloc.f90\
                         prec_const.f90 utils.f90
       interpolation.o : interpolation.f90 geometry.f90 MPIutils.f90 prec_const.f90\
                         spline1d.f90 utils.f90
        coord_system.o : coord_system.f90 geometry.f90 globals.f90 mem_alloc.f90\
                         prec_const.f90
         integration.o : integration.f90 coord_system.f90 geometry.f90 globals.f90\
                         OMPutils.f90 prec_const.f90 spline1d.f90
        read_profile.o : read_profile.f90 geometry.f90 globals.f90 prec_const.f90\
                         spline1d.f90 utils.f90
        init_profile.o : init_profile.f90 geometry.f90 globals.f90 mem_alloc.f90\
                         prec_const.f90 read_profile.f90 spline1d.f90 utils.f90
        init_current.o : init_current.f90 coord_system.f90 geometry.f90 globals.f90\
                         init_profile.f90 mem_alloc.f90 prec_const.f90
       init_magnetic.o : init_magnetic.f90 coord_system.f90 geometry.f90 globals.f90\
                         init_current.f90 init_profile.f90 mem_alloc.f90 prec_const.f90\
                         utils.f90
           fft_NRF90.o : fft_NRF90.f90 mem_alloc.f90 prec_const.f90
                 fft.o : fft.f90 fft_NRF90.f90 mem_alloc.f90 prec_const.f90
         gyroaverage.o : gyroaverage.f90 clock.f90 fft.f90 geometry.f90\
                         globals.f90 mem_alloc.f90 OMPutils.f90 prec_const.f90
            fequil4d.o : fequil4d.f90 geometry.f90 globals.f90 gyroaverage.f90\
                         init_magnetic.f90 init_profile.f90 mem_alloc.f90 MPIutils.f90\
                         prec_const.f90 spline1d.f90 utils.f90
         fdistribu5d.o : fdistribu5d.f90 bessel.f90 fequil4d.f90 geometry.f90\
                         globals.f90 init_magnetic.f90 init_profile.f90 local_spline.f90\
                         mem_alloc.f90 MPIutils.f90 prec_const.f90 spline1d.f90\
                         utils.f90
             poisson.o : poisson.f90 clock.f90 coord_system.f90 fdistribu5d.f90\
                         fft.f90 fft_NRF90.f90 geometry.f90 globals.f90\
                         gyroaverage.f90 init_magnetic.f90 init_profile.f90 integration.f90\
                         mem_alloc.f90 OMPutils.f90 prec_const.f90
              efield.o : efield.f90 geometry.f90 globals.f90 mem_alloc.f90\
                         OMPutils.f90 prec_const.f90 spline1d.f90
              filter.o : filter.f90 fdistribu5d.f90 fft.f90 globals.f90\
                         MPIutils.f90 OMPutils.f90 prec_const.f90
         advec1D_BSL.o : advec1D_BSL.f90 clock.f90 fdistribu5d.f90 geometry.f90\
                         globals.f90 init_magnetic.f90 init_profile.f90 integration.f90\
                         interpolation.f90 mem_alloc.f90 OMPutils.f90
         advec2D_BSL.o : advec2D_BSL.f90 clock.f90 fdistribu5d.f90 filter.f90\
                         geometry.f90 globals.f90 init_current.f90 init_magnetic.f90\
                         init_profile.f90 interpolation.f90 local_spline.f90 MPIutils.f90\
                         OMPutils.f90 prec_const.f90 spline1d.f90 utils.f90
              diag3D.o : diag3D.f90 clock.f90 f5D_saving.f90 fdistribu5d.f90\
                         globals.f90 gyroaverage.f90 init_current.f90 init_magnetic.f90\
                         integration.f90 OMPutils.f90 poisson.f90 prec_const.f90\
                         resu3D_saving.f90
      Pcross_section.o : Pcross_section.f90 fdistribu5d.f90 globals.f90 interpolation.f90\
                         mem_alloc.f90 MPIutils.f90 poisson.f90 prec_const.f90\
                         utils.f90
             physics.o : physics.f90 clock.f90 coord_system.f90 efield.f90\
                         fdistribu5d.f90 geometry.f90 globals.f90 gyroaverage.f90\
                         init_current.f90 init_magnetic.f90 init_profile.f90 mem_alloc.f90\
                         OMPutils.f90 prec_const.f90 utils.f90
          read_write.o : read_write.f90
       output_saving.o : output_saving.f90
         diagnostics.o : diagnostics.f90 clock.f90 coord_system.f90 diag3D.f90\
                         f5D_saving.f90 fdistribu5d.f90 geometry.f90 globals.f90\
                         gyroaverage.f90 init_current.f90 init_magnetic.f90 init_profile.f90\
                         integration.f90 mem_alloc.f90 MPIutils.f90 output_saving.f90\
                         Pcross_section.f90 physics.f90 poisson.f90 prec_const.f90\
                         read_write.f90 resu3D_saving.f90 utils.f90
      vlasov_solving.o : vlasov_solving.f90 advec1D_BSL.f90 advec2D_BSL.f90 clock.f90\
                         coord_system.f90 diagnostics.f90 efield.f90 fdistribu5d.f90\
                         filter.f90 geometry.f90 globals.f90 gyroaverage.f90\
                         init_current.f90 init_magnetic.f90 init_profile.f90 MPIutils.f90\
                         Pcross_section.f90 poisson.f90 prec_const.f90 read_write.f90
              gysela.o : gysela.f90 advec2D_BSL.f90 clock.f90 coord_system.f90\
                         diagnostics.f90 efield.f90 f5D_saving.f90 fdistribu5d.f90\
                         geometry.f90 globals.f90 gyroaverage.f90 gys_alloc.f90\
                         init_current.f90 init_magnetic.f90 init_profile.f90 integration.f90\
                         interpolation.f90 MPIutils.f90 OMPutils.f90 output_saving.f90\
                         Pcross_section.f90 physics.f90 poisson.f90 read_write.f90\
                         vlasov_solving.f90
$(OBJSOMP) :
	$(F90) -c $(CE_DFLAGS) $(OMPFLAG) $(R3FLAG) $(NOHDF5) $(TIFLAGS) $(MATRICIALFLAG) $(F90FLAGS) $(HDF5INCLUDE) $< -o $@

%.o: %.f90
	$(F90) -c $(CE_DFLAGS) $(R3FLAG) $(NOHDF5) $(TIFLAGS) $(MATRICIALFLAG) $(F90FLAGS) $(HDF5INCLUDE) $< -o $@

r3_ctlk.o:r3_ctlk.c
	$(CC) $(CFLAGS) -c r3_ctlk.c

ctime.o:ctime.c
	$(CC) $(CFLAGS) -c ctime.c


distclean: clean
	rm -f $(EXEC) $(EXEC_TI) make.dep

clean:	
	rm -f $(OBJS) *.o *~ *.mod  core
