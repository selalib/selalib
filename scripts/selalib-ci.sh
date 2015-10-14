#!/bin/bash
 
HOST=`hostname -s`
echo "${HOST}"
ARCH=`uname -s`
echo "ARCH:${ARCH}"
git log --date-order --date=short | \
sed -e '/^commit.*$/d' | \
awk '/^Author/ {sub(/\\$/,""); getline t; print $0 t; next}; 1' | \
sed -e 's/^Author: //g' | \
sed -e 's/>Date:   \([0-9]*-[0-9]*-[0-9]*\)/>\t\1/g' | \
sed -e 's/^\(.*\) \(\)\t\(.*\)/\3    \1    \2/g' > ChangeLog

if [[ `hostname` == "irma-hpc" ]]; then
  source /opt/intel/composerxe/bin/compilervars.sh intel64
  source /opt/intel/mkl/bin/mklvars.sh intel64
  source /opt/intel/impi/4.1.0.024/intel64/bin/mpivars.sh
  export FC=ifort
  export CC=icc
  export I_MPI_F90=ifort
  export I_MPI_CC=icc
  export HDF5_ROOT=/opt/local
  export FFTW_ROOT=/opt/local
fi

rm -rf build
mkdir build
cd build; {
${CMAKE} \
	-DCMAKE_BUILD_TYPE=Release \
	-DHDF5_PARALLEL_ENABLED=ON \
	${WORKSPACE} 
make Experimental
}; cd -

rm -rf build

exit 0


mkdir build
cd build; {
cmake ../
make NightlyUpdate
make NightlyConfigure
make NightlyBuild
make NightlyTest || true
if [ -f Testing/TAG ] ; then
   xsltproc ../cmake/ctest2junix.xsl Testing/`head -n 1 < Testing/TAG`/Test.xml > CTestResults.xml
fi
make NightlySubmit
}; cd -
rm -rf build
exit 0
