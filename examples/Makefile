SLL_ROOT = /Users/navaro/Codes/selalib/usr
F90 = gfortran
OPT = -O3 #-xHost -ip -openmp
LD  = gfortran
   
F90FLAGS = -I${SLL_ROOT}/include/selalib
LDFLAGS = -L${SLL_ROOT}/lib 
LIBS = -lselalib -ldfftpack

all: $(PROG)

landau.exe: landau.F90
	$(LD) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

run:
	./landau | gnuplot
 
clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .F90

.F90.o:
	$(F90) $(F90FLAGS) -c $<

.mod.o:
	$(F90) $(F90FLAGS) -c $*.F90
