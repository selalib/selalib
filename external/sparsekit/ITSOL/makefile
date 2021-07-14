
F77 = gfortran

#F77 = gggg
FFLAGS = -g -Wall 

LIBS =  ../libskit.a ../UNSUPP/BLAS1/blas1.o

riters.ex:  riters.o iters.o ilut.o itaux.o $(LIBS)
	$(F77) -o riters.ex riters.o itaux.o $(LIBS)

rilut.ex: rilut.o ilut.o iters.o itaux.o $(LIBS)
	$(F77) -o rilut.ex rilut.o itaux.o $(LIBS)

riter2.ex: riter2.o iters.o ilut.o itaux.o $(LIBS)
	$(F77) -o riter2.ex riter2.o itaux.o $(LIBS)

%.o : %.f
	$(F77) $(FFLAGS) -c $<


clean:
	rm -f *.o *.ex core *.trace


../UNSUPP/BLAS1/blas1.o:
	(cd ../UNSUPP/BLAS1; $(F77) $(FFLAGS) -c blas1.f)

../libskit.a:
	(cd ..; $(MAKE) $(MAKEFLAGS) libskit.a)
