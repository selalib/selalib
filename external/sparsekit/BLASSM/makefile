FFLAGS = -DUNIX -O3 -g 
F77 = gfortran
LINK = -L ../. -lskit

#F77 = cf77
#FFLAGS = -Wf "-dp"

mvec.ex: rmatvec.o ../MATGEN/FDIF/functns.o ../libskit.a
	$(F77) $(FFLAGS) -o mvec.ex rmatvec.o ../MATGEN/FDIF/functns.o $(LINK)

tester.ex: tester.o ../MATGEN/FDIF/functns.o 
	$(F77) $(FFLAGS) -o tester.ex tester.o ../MATGEN/FDIF/functns.o $(LINK)
clean:
	rm -f *.o *.ex core *.trace fort.* ftn??

../MATGEN/FDIF/functns.o:
	(cd ../MATGEN/FDIF; $(F77) $(FFLAGS) -c functns.f)

../libskit.a:
	(cd ..; $(MAKE) $(MAKEFLAGS) libskit.a)

%.o : %.f
	$(F77) $(FFLAGS) -o $@ -c $<

