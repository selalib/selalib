#!/bin/bash
 
# Every night, each entry following '#PIPOL' is passed as argument to
# a pipol-sub command:

# There is two ways to make a reservation via pipol-sub:
# --> Add a line per image with the following syntaxe:
# '#PIPOL esn amd64-linux-debian-testing.dd none 2:00 --silent --user'
#  (1)    (2) (3)                           (4)  (5)  (6)      (7)
     
# (1) comment token known by pipol
# (2) a reservation as soon as possible (at night, reservations are starting at 20h30)
# (3) the system image
# (4) none as we do not want a particular host
# (5) the reservation is for two hours
# (6) do not send reservation emails
# (7) run as user (You can only run the script as root for the Windows Images).

# Use entries only if you need them!

# --> The second way is to use pipol-sub shortcut:
# '#PIPOL debian 2:00 --root`
#  (1)         (2)       (3)    (4)

# (1)  Comment token known by pipol
# (2)  The string to matched with. 
#	   In this example the line will reserve and deploy all the Debian images.
#        For example, the line: #PIPOL linux 2:00 --root will reserve and deploy all the linux images.
#        For more informations, please read the documentation: http://pipol.inrialpes.fr/documentation/manual/user-manual/user-manual.html#ex-cutions-en-batch-mode-simplifi
# (3)  The reservation is for two hours
# (4)  Some options, please read the documentation for more details

#PIPOL amd64_2010-linux-debian-testing 1:00 
#PIPOL amd64-linux-debian-testing 1:00 
#PIPOL amd64_2010-linux-fedora-core16 1:00
#PIPOL amd64_2010-linux-fedora-core17 1:00
#PIPOL amd64-linux-fedora-core16 1:00
#PIPOL amd64-linux-ubuntu-oneiric 1:00
#PIPOL amd64-linux-ubuntu-precise 1:00
#PIPOL amd64-linux-ubuntu-quantal 1:00
#PIPOL x86_64_mac-mac-osx-server-snow-leopard 2:00

HOST=`hostname -s`
echo "${HOST}"
ARCH=`uname -s`

if [[ $PIPOL_IMAGE_NAME ]]; then
   echo "PIPOL_IMAGE_NAME:${PIPOL_IMAGE_NAME}"
else
   echo "ARCH:${ARCH}"
fi

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

if [[ $PIPOL_WDIR ]] ; then
   echo "PIPOL_WDIR:$PIPOL_WDIR"
   WORKDIR=$PIPOL_WDIR
   echo "PIPOL_HOMEDIR:$PIPOL_HOMEDIR"
   echo "PIPOL_IMAGE_NAME:${PIPOL_IMAGE_NAME}"
   HOMEDIR=$PIPOL_HOMEDIR
elif [[ "$HOST" = *-navaro* ]]; then
   WORKDIR="/tmp"
   HOMEDIR="${HOME}/Codes"
elif [[ "$HOST" =~ "vpn-irma" ]]; then
   HOMEDIR="${HOME}/Codes"
   WORKDIR="/tmp"
elif [[ "$HOST" =~ "doct" || "$HOST" =~ "reserve" ]]; then
   WORKDIR="/Users/irma"
   HOMEDIR="/Users/irma"
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

if [[ $PIPOL_IMAGE_NAME == *-fedora-* ]]; then
   echo "Fedora core"
   echo "Try to load module "
   source /etc/profile.d/modules.sh 
   module load openmpi-x86_64
fi

if [ "$PIPOL_IMAGE_NAME" = "x86_64_mac-mac-osx-server-snow-leopard" ]; then
   echo "Configuration pour MACOS SNOWLEOPARD"
   export PATH=/usr/local/bin:${PATH}
   export CC=/usr/local/bin/gcc
   export CXX=/usr/local/bin/g++
   export FC=/usr/local/bin/gfortran
   export OMPI_FC=/usr/local/bin/gfortran
   export OMPI_CC=/usr/local/bin/gcc
   export HDF5_ROOT=/usr/local
   export FFTW_ROOT=/usr/local
fi

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
  alias cmake='/usr/local/bin/cmake'
fi

mkdir build
cd build; {
cmake ${HOMEDIR}/selalib -DCMAKE_BUILD_TYPE=Release
make NightlyUpdate
make NightlyConfigure
make NightlyBuild
make NightlyTest
make NightlySubmit
}; cd -

rm -rf build

exit 0
