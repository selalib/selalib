#!/bin/bash -e

# ---
# Script to build selalib on MPCDF systems, in particular on COBRA and DRACO.
#
# Important: This script needs to be run *in place* from the ./scripts/ folder
# in order to be able to correctly determine the source location.  The build is
# performed in a separate directory ($BUILD_DIR, see below) which defaults to
# ~/selalib_obj.
#
# Note that various switches and flags can be set at the script invocation, e.g. via
#   $ DO_CONF=0 JMAKE=16 ./compile_mpcdf.sh
# without the need to edit this script.
# ---


# run the cmake configure step
DO_CONF=${DO_CONF:=1}

# run the actual build
DO_BUILD=${DO_BUILD:=1}

# start with an empty object directory
DO_CLEAN=${DO_CLEAN:=0}

# define the CMAKE build type
CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:="Release"}

# build the external package "fmempool"
USE_FMEMPOOL=${USE_FMEMPOOL:=OFF}
#USE_FMEMPOOL=${USE_FMEMPOOL:=ON}

# number of processors to be used for a parallel build
JMAKE=${JMAKE:=8}

# add custom preprocessor macros (e.g. toggle 32 bit halos) as follows
#MACROS="-DUSE_HALO_REAL32"
#MACROS="-DDISABLE_CACHE_BLOCKING"
MACROS=${MACROS:=""}

# BUILD_DIR_PRESET=$HOME/selalib_obj
# With many small files /tmp gives much better performance compared to the
# home directory which often resides on a large shared GPFS file system.
#umask 0077
BUILD_DIR_PRESET=/tmp/$USER/selalib_obj
#export ZFP_ROOT=/afs/ipp/home/k/khr/soft/amd64_sles11/opt/zfp/0.5.1-intel-16
BUILD_DIR=${BUILD_DIR:="$BUILD_DIR_PRESET"}

# ifort: add the "-ipo-separate" flag (interprocedural optimization)
# USE_IPO=${USE_IPO:=0}
USE_IPO=${USE_IPO:=1}

# ifort: add the "-fp-model precise" flag (better reproducibility, but deteriorates vectorization and speed)
USE_FP_PRECISE=${USE_FP_PRECISE:=0}

# ifort: align arrays to 64 byte boundaries to improve vectorization
USE_ALIGNMENT=${USE_ALIGNMENT:=1}

# ifort: add various error checking compiler flags (slowdown!)
USE_CHECKS=${USE_CHECKS:=0}

USE_CLAPP=${USE_CLAPP:=0}
if [ x"$USE_CLAPP" != x"0" ]; then
    export CLAPP_DIR=${CLAPP_DIR:=$HOME/clapp/usr}
fi

# --- end of configuration section ---


# automatically determine the absolute location of the selalib source tree
SOURCE_DIRECTORY=`pwd`/../
SOURCE_DIRECTORY=`readlink -f $SOURCE_DIRECTORY`
LOGFILE="$BUILD_DIR/build.log"
SEP="-----------------------------------------------------------------------------"


if [ ! -d "../scripts" ]
then
    echo $SEP
    echo "Error: This script must be run from the './scripts' subdirectory of SeLaLib,"
    echo "i.e. typically as follows: './compile_mpcdf.sh'"
    echo $SEP
    exit 1
fi


if [ x"$DO_CLEAN" == x"1" ] && [ -e "$BUILD_DIR" ]
then
    SUFFIX=`date +%Y%m%d%H%M%S`
    BACKUP="${BUILD_DIR}.${SUFFIX}"
    echo $SEP
    echo "Moving existing build directory to ${BACKUP}"
    echo $SEP
    mv "$BUILD_DIR" "$BACKUP"
fi


mkdir -p $BUILD_DIR
# apply redirection to log file
{
    # --- modules and machine-dependent optimization flags
    source /etc/profile.d/modules.sh
    module purge
    # --- We're currently facing compiler errors with the 2019 version of Intel:
    # module load intel
    # module load impi
    # module load mkl
    module load intel/18.0.5
    module load impi/2018.4
    module load mkl/2018.4
    # ---
    module load fftw-serial
    module load hdf5-mpi
    module load cmake
    module load anaconda/3
    module load git

    INTEL_OPT_FLAGS="-O3 -xHost"
    # Skylake machines
    if [ x"$CLUSTER" == x"COBRA" ] || [ x"$CLUSTER" == x"SAKURA" ] || [ x"$CLUSTER" == x"TALOS" ]; then
        # high usage of AVX512 registers seems slightly beneficial
        INTEL_OPT_FLAGS="$INTEL_OPT_FLAGS -qopt-zmm-usage=high"
    fi
    INTEL_OPT_FLAGS="$INTEL_OPT_FLAGS -nowarn"

    # align to 64 byte boundaries to enable vectorization
    if [ x"$USE_ALIGNMENT" == x"1" ]
    then
        INTEL_F90_FLAGS="$INTEL_F90_FLAGS -align array64byte"
    fi

    # WARNING : To get accurate results (unit tests)
    # we need to turn off aggressive FP optimizations as follows:
    # (fp-model precise costs about 2%-5%)
    if [ x"$USE_FP_PRECISE" == x"1" ]
    then
        INTEL_OPT_FLAGS="$INTEL_OPT_FLAGS -fp-model precise"
    fi

    # Inter-file interprocedural optimization accelerates by nearly 10%,
    # but makes compilation time significantly longer and sometimes produces internal compiler errors (seen with Intel 2017)
    if [ x"$USE_IPO" == x"1" ]
    then
        INTEL_OPT_FLAGS="$INTEL_OPT_FLAGS -ipo-separate"
    fi

    if [ x"$USE_CHECKS" == x"1" ]
    then
        INTEL_OPT_FLAGS="$INTEL_OPT_FLAGS -check"
        INTEL_OPT_FLAGS="$INTEL_OPT_FLAGS -fpe0 -traceback"
    fi

    # insert rarely-needed special flags such as "-qopt-report" manually
    if [ x"$INTEL_SPECIAL_FLAGS" != "" ]
    then
        INTEL_OPT_FLAGS="${INTEL_OPT_FLAGS} ${INTEL_SPECIAL_FLAGS}"
    fi

    # boundary checking does only work for Fortran arrays (icc would complain)
    # INTEL_F90_FLAGS="-check bounds"

    echo $SEP
    echo "SeLaLib build configuration"
    echo $SEP
    echo "DO_CONF=${DO_CONF}"
    echo "DO_BUILD=${DO_BUILD}"
    echo "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
    echo "USE_FMEMPOOL=${USE_FMEMPOOL}"
    echo "JMAKE=${JMAKE}"
    echo "MACROS=${MACROS}"
    echo "INTEL_OPT_FLAGS=${INTEL_OPT_FLAGS}"
    echo "INTEL_F90_FLAGS=${INTEL_F90_FLAGS}"
    echo "BUILD_DIR=${BUILD_DIR}"
    echo "SOURCE_DIRECTORY=${SOURCE_DIRECTORY}"
    echo "LOGFILE=${LOGFILE}"
    echo "CLUSTER=${CLUSTER}"
    echo $SEP
    echo "Proceeding to build ..."
    echo

    module list
    echo

    # --- set enviroment variables to actually use the Intel toolchain
    export FC=ifort
    export CC=icc
    export CXX=icpc

    mkdir -p $BUILD_DIR
    cd $BUILD_DIR

    if [ x"$DO_CONF" == x"1" ]
    then
      cmake $SOURCE_DIRECTORY \
        -DUSE_FMEMPOOL:BOOL=${USE_FMEMPOOL} \
        -DOPENMP_ENABLED:BOOL=ON \
        -DUSE_MKL:BOOL=ON \
        -DHDF5_PARALLEL_ENABLED:BOOL=ON \
        -DCMAKE_BUILD_TYPE:STRING="${CMAKE_BUILD_TYPE}" \
        -DMPI_C_COMPILER:STRING="mpiicc" \
        -DMPI_CXX_COMPILER:STRING="mpiicpc" \
        -DMPI_Fortran_COMPILER:STRING="mpiifort" \
        -DFORCE_Fortran_FLAGS_RELEASE:STRING="$INTEL_OPT_FLAGS" \
        -DCMAKE_C_FLAGS:STRING="-g $MACROS $INTEL_OPT_FLAGS" \
        -DCMAKE_CXX_FLAGS:STRING="-g $MACROS $INTEL_OPT_FLAGS" \
        -DCMAKE_Fortran_FLAGS:STRING="-g $INTEL_F90_FLAGS $MACROS" \
        -DCLAPP:BOOL=${USE_CLAPP}
    else
      echo "skipping configure step ..."
    fi

    if [ x"$DO_BUILD" == x"1" ]
    then
      gmake -j${JMAKE} VERBOSE=1
    else
      echo "skipping make step ..."
    fi

    echo $SEP
    echo "OK!"
    echo $SEP
} 2>&1 | tee "$LOGFILE"
