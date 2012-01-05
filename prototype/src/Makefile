
default:
	(cd memory; scons)
	(cd precision; scons)
	(cd assert; scons)
	(cd constants; scons)
	(cd utilities; scons)
	(cd diagnostics; scons)
	(cd non_linear_solvers; scons)
	(cd splines; scons)
        (cd fftpack; make)
	(cd fft; scons)
	(cd collective; scons CPATH='.' LPATH='.')
	(cd remap; scons)
	(cd timer; scons)
	(cd ode_solvers; scons)
	(cd mesh_types; scons)
	(cd distribution_function; scons)

install_fftw:
	(cd fftw; make)

install_fftpack:
	(cd fftpack; make)

clean:
	(cd constants; scons -c)
	(cd precision; scons -c)	
	(cd memory; scons -c)
	(cd assert; scons -c)
	(cd non_linear_solvers; scons -c)
	(cd splines; scons -c)
	(cd utilities; scons -c)
	(cd fft; scons -c)
	(cd diagnostics; scons -c)
	(cd integration; scons -c)
	(cd remap; scons -c)
	(cd collective; scons -c)
	(cd timer; scons -c)
	(cd mesh_types; scons -c)
	(cd ode_solvers; scons -c)
	(cd distribution_function; scons -c)
