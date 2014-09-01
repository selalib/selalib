module purge
module load cmake
module load gnu-env
module load openmpi/1.6.3_gnu47
module load fftw/3.3.3_gnu47     
module load hdf5/1.8.10_gnu47_openmpi  
module load lapack
export HDF5_ROOT=${HDF5_ROOT_DIR}
export FFTW_ROOT=${FFTW_INC_DIR}/../
module show openmpi/1.6.3_gnu47
cmake -DHDF5_PARALLEL_ENABLED=ON SELALIB_DIRECTORY/prototype/src
