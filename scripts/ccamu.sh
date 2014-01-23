module purge
module load cmake/2.8.7 hdf5/gcc47/1.8.8 fftw/gcc47/3.3
export HDF5_ROOT=/softs/hdf5/gcc47/1.8.8/
export FFTW_ROOT=/softs/fftw/gcc47/3.3/
cd /scratch/pnavaro
mkdir selalib-build
cd selalib-build; {
cmake ${HOME}/selalib -DCMAKE_BUILD_TYPE=Release -DOPTIONS_FILE=${HOME}/selalib/cmake/ccamu.cmake
make NightlyUpdate
make NightlyConfigure
make NightlyBuild
/usr/bin/oarsub -S ${HOME}/selalib/scripts/ccamu-gcc.oar
make NightlySubmit
}; cd -
rm -rf /scratch/pnavaro/selalib-build
module purge
module load hdf5/intel/1.8.8-cemracs
module load fftw/intel/3.3
module load intel/13.1.2
module load cmake/2.8.7
export HDF5_ROOT=/softs/cemracs/selalib/intel
export FFTW_ROOT=/softs/fftw/intel/3.3
cd /scratch/pnavaro
mkdir selalib-intel
cd selalib-intel; {
cmake ${HOME}/selalib -DOPTIONS_FILE=${HOME}/selalib/cmake/ccamu-intel.cmake
make NightlyUpdate
make NightlyConfigure
make NightlyBuild
/usr/bin/oarsub -S ${HOME}/selalib/scripts/ccamu-intel.oar
make NightlySubmit
}; cd -
rm -rf /scratch/pnavaro/selalib-intel
