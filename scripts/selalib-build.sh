#!/bin/bash
 
HOST=`hostname -s`
echo "${HOST}"
ARCH=`uname -s`
CMAKE='cmake'
echo "ARCH:${ARCH}"

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

if [[ "$HOST" = *-navaro* ]]; then
   WORKDIR="/tmp"
   HOMEDIR="${HOME}/Codes"
elif [[ "$HOST" =~ "vpn-irma" ]]; then
   HOMEDIR="${HOME}/Codes"
   WORKDIR="/tmp"
elif [[ "$HOST" =~ "doct" || "$HOST" =~ "reserve" ]]; then
   WORKDIR="/Users/irma"
   HOMEDIR="/Users/irma"
elif [[ "$HOST" = irma-gpu3 ]]; then
   WORKDIR="/scratch/navaro"
   HOMEDIR=${HOME}
elif [[ "$HOST" = *irma-* ]]; then
   WORKDIR="/scratch/navaro"
   HOMEDIR=${HOME}/Vlasov
else
   HOMEDIR="./"
   WORKDIR="./"
fi

echo "HOMEDIR:$HOMEDIR"
echo "WORDIR:$WORKDIR"

if [[ $(($(stat -f --format="%a*%S" $WORKDIR))) == 0 ]]; then
   cd /tmp
else
   cd $WORKDIR
fi

cd ${HOMEDIR}/selalib; {
   git log --date-order --date=short | \
   sed -e '/^commit.*$/d' | \
   awk '/^Author/ {sub(/\\$/,""); getline t; print $0 t; next}; 1' | \
   sed -e 's/^Author: //g' | \
   sed -e 's/>Date:   \([0-9]*-[0-9]*-[0-9]*\)/>\t\1/g' | \
   sed -e 's/^\(.*\) \(\)\t\(.*\)/\3    \1    \2/g' > ChangeLog
}; cd -

if [[ `hostname -s` == "irma-4600" ]]; then
  source /opt/intel/composerxe/bin/compilervars.sh intel64
  source /opt/intel/mkl/bin/mklvars.sh intel64
  source /opt/intel/impi/4.1.1.036/bin64/mpivars.sh 
  export CC=icc
  export I_MPI_F90=ifort
  export I_MPI_CC=icc
  export HDF5_ROOT=/opt/local
  export FFTW_ROOT=/opt/local
  export MKLROOT=/opt/intel/composer_xe_2013.4.183/mkl
  export CMAKE=/usr/local/bin/cmake
fi

mkdir build
cd build; {
${CMAKE} \
	-DCMAKE_BUILD_TYPE=Release \
	-DHDF5_PARALLEL_ENABLED=ON \
	${HOMEDIR}/selalib/prototype/src 
make NightlyUpdate
make NightlyConfigure
make NightlyBuild
make NightlyTest
make NightlySubmit
}; cd -

rm -rf build

exit 0
