
default:
	( gfortran -c dfftpack.f )
	( gfortran -c fftpack.f90 )
	( gfortran -c fftpack_module.F90 )
	( ar -r libfftpack.a fftpack_module.o fftpack.o dfftpack.o )
	#( cd ../fft ; bash install.sh 1 >/dev/null)

clean:
	( rm -f *.o )
	( rm -f *.a )
	( rm -f *.mod )
	#( cd ../fft ; bash install.sh 1 >/dev/null)
