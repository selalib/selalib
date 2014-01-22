#! /bin/bash
#SBATCH -p public
#SBATCH -A math
#SBATCH -n 16     
#SBATCH -N 2-4     
#SBATCH -t 01:00:00 
#SBATCH --mem=1024   
#SBATCH --mail-type=END 
#SBATCH --mail-user=navaro@unistra.fr

module purge
module load batch/slurm
module load compilers/intel13
module load mpi/openmpi-1.6.i13.threaded
module load libs/mkl13
module load libs/hdf5-1.8
export HDF5_ROOT=/usr/local/hdf5-1.8.5
export FFTW_ROOT=/usr/local/fftw-3.3.2
source /opt/intel/composer_xe_2013.2.146/mkl/bin/mklvars.sh intel64
cd /workdir/math/navaro/build/
make test
