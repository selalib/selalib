source ${HOME}/.bashrc
export HDF5_ROOT=/usr/local/hdf5-1.8.5
export FFTW_ROOT=/usr/local/fftw-3.3.2
source /opt/intel/composer_xe_2013.2.146/mkl/bin/mklvars.sh intel64
cd /workdir/math/navaro/all
make NightlyUpdate
make NightlyConfigure
make NightlyBuild
/usr/local/slurm/bin/sbatch ${HOME}/selalib/scripts/hpc-login.sl
make NightlySubmit
