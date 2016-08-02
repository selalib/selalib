module purge
module load intel-env/15.0.0 
module load openmpi/1.6.3_intel15.0.0 
module load hdf5/1.8.16_intel_openmpi
module load cmake/3.4.3 
module load mkl/11.2 
module load python/epd
export HDF5_ROOT=${HDF5_ROOT_DIR}
export FFTW_ROOT=${FFTW_ROOT_DIR}/../
cmake ~/selalib -DHDF5_PARALLEL_ENABLED=ON -DUSE_MKL=ON -DCMAKE_Fortran_COMPILER=h5pfc
