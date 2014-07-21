source /etc/profile
module load cmake
module load hdf5/mpi
module load fftw 
module load lapack
export HDF5_ROOT=/bglocal/cn/pub/HDF5/1.8.9/par
export FFTW_ROOT=/bglocal/cn/pub/FFTW/3.3.3/
export SELALIB_DIR=${HOME}/selalib
mkdir ${workdir}/build
cd ${workdir}/build
cmake -DCMAKE_TOOLCHAIN_FILE=${SELALIB_DIR}/cmake/turing-toolchain.cmake \
      -DOPTIONS_FILE=${SELALIB_DIR}/cmake/turing.cmake \
      ${SELALIB_DIR}/prototype/src
