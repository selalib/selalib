
default:
	( tar zxvf fftw-3.3.tar.gz )
	#( export TOTO=`pwd` )
	#( TOTO=${TOTO}/ )
	( ./fftw-3.3/configure --prefix=`pwd` ; cd fftw-3.3/ ; make ; make install )
	#( gfortran -c -L./fftw-3.3/lib -I./fftw-3.3/include fftw.f03 -lfftw3 -lm )
	( gfortran -c -L./lib -I./include fftw.f03 )
	#( ar -r libFFTW3.a fftw.o )
	#( cd ../fft ; bash install.sh 1 >/dev/null)

clean:
	( cd fftw-3.3/ ; make uninstall )
	( rm fftw.o ; rm fftw3.mod )
	#( cd ../fft ; bash install.sh 1 >/dev/null)