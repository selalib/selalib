module load cmake
module load intel-env
module load mkl
module load intelmpi
module load fftw
module load hdf5
export HDF5_ROOT=${HDF5_ROOT_DIR}
export FFTW_ROOT=${FFTW_ROOT_DIR}/../
cmake -DCMAKE_Fortran_COMPILER=ifort \
      -DUSE_MKL=ON \
      -DHDF5_PARALLEL_ENABLED=ON SELALIB_DIRECTORY/prototype/src
