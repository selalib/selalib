#!/bin/bash
#SBATCH -J SELALIB                             # job name
#SBATCH -A selavlas                            # account name
#SBATCH -N 512                                 # number of nodes
#SBATCH -n 8192                                # number of processes
###SBATCH --place                              # make sure --place is turned off
#SBATCH -o %j.out                              # stdout file name (%j: job ID)
#SBATCH -e %j.err                              # stderr file name (%j: job ID)
#SBATCH -t 01:00:00                            # max run time (hh:mm:ss)

. /etc/profile.d/00-modules.sh
module purge
#module load cmake intel
module load intelmpi hdf5_p fftw srun
export FC=ifort
export CC=icc
export CXX=icpc
export ARCH=helios


#######################################
JOBDIR=$SLURM_JOB_ID
#This assumes you have selalib in your home path
#and a spearate build directory
SELALIB_BUILD=$HOME/selalib_build
SELALIB=$HOME/selalib
#Binary to be used
SIMULATION_BINARY=test_general_pif
SIMULATION_INPUTFILE=$SELALIB/src/simulations/simulations_parallel/sim_general_vlasov_poisson_pif/landau_magnetized.nml
######################################

#Create job directory in Work directory
echo "Working directory is: $WORKDIR"
cd $WORKDIR
mkdir $JOBDIR
cd $JOBDIR

#Copy binary into jobdir
cp $SELALIB_BUILD/bin/$SELALIB_BINARY ./

srun $SELALIB_BINARY $SIMULATION_INPUTFILE
