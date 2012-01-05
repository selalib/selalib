
TOTO=$(shell pwd)

default:
	( tar -xzf fftw-3.3.tar.gz )
	( cd fftw-3.3/ ; ./configure --prefix=TOTO ; make ; make install )
	( gfortran -c -L./lib -I./include fftw.f03 )

clean:
	( cd fftw-3.3/ ; make uninstall )
	( rm -f fftw.o ; rm -f fftw3.mod )
