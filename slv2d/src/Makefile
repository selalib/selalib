.SUFFIXES: .o .f90 .F90 .f .swp
#
#  Makefile for vp2d on Solaris
#
OBJECTS =  precision.o module_mpi.o geometry.o diagnostiquesp.o	\
	   fft_module.o poisson2ddir.o poisson2dpp.o \
	   splinepx.o splinepy.o splinepp.o splinenn.o \
 	   vlasov2d.o vlasov1d.o transpose.o maxwell2dfdtd.o diagnostiquesm.o \

ARCHIVES = libdsine.a

## irma-hpc ##
F90 = mpif90 -c -cpp -fbounds-check
LD  = mpif90
AR  = ar -r
HDF5 = /usr/local#silo est au meme endroit
INCL = -I$(HDF5)/include 
OPT = -m64
F90FLAGS = $(OPT) $(INCL) -D_MPI -D_FFTPACK
LIB = -llapack -L$(HDF5)/lib -lhdf5 -lsiloh5 -ldfftpack -lz -lgfortran -lmpi -lsupc++
#LIB = -xlic_lib=sunperf -L$(HDF5)/lib -lhdf5 -lsiloh5 -R/usr/local/hdf5/lib

## gateway EFDA ##
#F90 = mpif90 -c
#LD  = mpif90 
#AR  = ar -r
#HDF5=/afs/efda-itm.eu/project/switm/visit/hdf5/1.8.0/linux
#SILO=/afs/efda-itm.eu/project/switm/visit/silo
#INCL = -I$(HDF5)/include -I$(SILO)/include 
#F90FLAGS = $(OPT) $(INCL)
#LIB = -L$(HDF5)/lib -lhdf5 -L$(SILO)/libsilo.a \
#      -L/afs/efda-itm.eu/project/switm/libAMD64 -lrfftw

## hpc CECPV ##
#F77 = ifort -c 
#F90 = mpif90 -c 
#LD  = mpif90 
#AR = ar -r
#FFLAGS = -O3 #-fast  
#F90FLAGS = -O3 #-fast 
#LIB = ../dfftpack/libdfftpack.a -lblas -llapack

.f90.o:	
	$(F90) $(F90FLAGS) $<

.mod.o:	
	$(F90) $(F90FLAGS) $*.f90

.F90.o:
	$(F90) $(F90FLAGS) $<  

.f.o:	
	$(F90) $(F90FLAGS) $<

default: vp2d vm2d

vm2d:	$(OBJECTS) vm2dinit.o vm2d.o $(ARCHIVES)
	$(LD) $(F90FLAGS) -o ../vm2d $^ $(LIB) 

vm2dinit.o:vm2dinit.f90

vp2d:	$(OBJECTS) vp2dinit.o vp2d.o $(ARCHIVES)
	$(LD) $(F90FLAGS) -o ../vp2d $^ $(LIB) 

vp2dinit.o:vp2dinit.F90
	$(F90) $(F90FLAGS) -D_MPI -D_FFTPACK $(INCL) vp2dinit.F90

mabo2d:	$(OBJECTS) mabo2dinit.o mabo2d.o $(ARCHIVES)
	$(LD) $(F90FLAGS) -o ../mabo2d $(OBJECTS) mabo2dinit.o mabo2d.o $(ARCHIVES) $(LIB)
 
libdsine.a: vsint.o vsinti.o 
	$(AR) libdsine.a vsint.o vsinti.o 

testp:  $(OBJECTS) test_electro2d.o
	$(LD) $(F90FLAGS) -o testp $(OBJECTS) test_electro2d.o $(ARCHIVES) $(LIB) 

clean:	
	rm -f *.o *~ *.mod  *.l *.T testp mabo2d vp2d vm2d

debug:	
	make F90FLAGS="-g -Rb -C" FFLAGS="-g -C"

perf:
	make F90FLAGS="-O3,unroll2,pipeline2 -eA -lapp"
