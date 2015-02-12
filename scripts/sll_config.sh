#!/bin/bash
EXTRA_ARGS=$@
ARCH="`uname -p`"
OS="`uname -s`"
HOSTNAME=`hostname -s`

WD="`pwd`"
# -------------------------------------------------------
# some useful pathnames
# change these if you want; leave everything else alone
# make sure path names are absolute
SLL_DIR="${HOME}/Codes/selalib"
SLL_BUILD="$WD/build_with_macports"
# -------------------------------------------------------

INFO( ) {
   echo "+---------------------------------------------"
   echo "+ $*"
   echo "+---------------------------------------------"
   [ $OS = macosx ] && say "$*"
}

DIE( ) { echo "***** ERROR ***** $*"; exit 1; }

USAGE( ) {
     echo "Usage:  ./sll_config [OPTIONS]"
     echo ''
     echo "    Installer for SeLaLib on MacOSX/Linux."
     echo ''
     echo "OPTIONS:"
     echo " -d SOURCE_PATH    Filesystem path to selalib directory"
     echo " -b BUILD_PATH     Filesystem path to build directory"
     echo " -h                Display this help and exit"
     echo " -v                Output version information and exit"
     echo ''
     echo "    Report bugs to <selalib-users@lists.gforge.inria.fr>"
  }

INFO "Operating System is $OS"
# make sure the OS is supported by this script
case "$OS" in
   Darwin)  OS="macosx"  ;;
   Linux*)  OS="linux"   ;;
        *)  DIE "unkown OS - you need to modify the installer script";;
esac

case $ARCH in
   i?86)  ARCH="i386";;
   x86_64)  ARCH="x86_64";;
        *)  uname -a | grep -q 'i.86'   && ARCH="i386"
            uname -a | grep -q 'x86_64' && ARCH="x86_64";;
esac


# get the command options, if any
while getopts ":d:b:vh" opt; do
   case $opt  in
      d ) SLL_DIR="$OPTARG";;
      b ) SLL_BUILD="$OPTARG";;
      v ) echo "SELALIB version $VERSION"; exit 0;;
      h ) USAGE; exit 0;;
     \? ) USAGE; DIE "bad option" ;;
   esac
done

INFO "Your hostname is : $HOSTNAME"
INFO "The selalib source directory is : $SLL_DIR"

if [ -d "$SLL_BUILD" ] ; then
   INFO  "$SLL_BUILD already exists - moving it to ${SLL_BUILD}_old"
   mv    "$SLL_BUILD" "${SLL_BUILD}_old"
fi
mkdir -p "$SLL_BUILD" ||  DIE "can-t create $SLL_BUILD"

INFO "make sure MAKE and C MAKE is installed"

case "$OS" in
   macosx)  CMAKE="/opt/local/bin/cmake"; MAKE="/usr/bin/make" ;;
   linux)   CMAKE="`which cmake`";   MAKE="`which make`" ;;
esac

[ -x "$CMAKE"  ] || DIE "missing executable file = $CMAKE"
[ -x "$MAKE" ] || DIE "missing executable file = $MAKE"

INFO "check OS & cpu type"

case "$ARCH" in
       i386) : ;;
     x86_64) : ;;
          *) echo " ***** unknown cpu type: $ARCH";;
esac
echo "...OS:     $OS"
echo "...ARCH:   $ARCH"

INFO "Change directory to $SLL_BUILD"
cd $SLL_BUILD

case "$OS" in
   macosx) 
      PORT="/opt/local/bin/port" 
      if [ -s $PORT ] ; then 
         INFO "macports is installed - configure the library"
         ${CMAKE} -DOPTIONS_FILE:FILEPATH=$SLL_DIR/cmake/macports_config.cmake \
             $EXTRA_ARGS ${SLL_DIR}/prototype/src 
      else 
         ${CMAKE} $EXTRA_ARGS ${SLL_DIR}/prototype/src 
      fi;;
   linux)  ${CMAKE} $EXTRA_ARGS ${SLL_DIR}/prototype/src ;;
esac

