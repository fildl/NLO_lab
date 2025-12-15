#!/bin/bash
#---------------------------------------------------------------------------- #
#   University    |   DIFA - Dept of Physics and Astrophysics 
#       of        |   Open Physics Hub
#    Bologna      |   (https://site.unibo.it/openphysicshub/en)
#----------------------------------------------------------------------------
#
# Usage
#   run a job:         sbatch run.sh
#   check processes:   slurmtop
#   delete a job:      scancel <jobID>   
#
#SBATCH --nodes=1                        ## number of nodes to be allocated
#SBATCH --tasks-per-node=2               ## number of tasks per node
#SBATCH --ntasks=2                       ## total number of tasks (nodes x tasks-per-node)
#SBATCH --mem-per-cpu=4G                 ## ram per cpu (to be tuned)
#SBATCH --job-name="filo-tides"          ## job name in the scheduler
#SBATCH --output=nemo_script.out         ## log file 
#SBATCH --error=nemo_script.err          ## err file

# --------------------------------------------------------------------------- #
# Modules setup and applications run
# --------------------------------------------------------------------------- #
module load ucx/1.13.1
module load mpi/openmpi/4.1.4 
module load terra/NEMO_env/1.0.1

# --------------------------------------------------------------------------- #
# Create a worning dir (clean)
# orig_path is the EXP00 for the specific experiment
# i.e. 
# orig_path='/home/PERSONALE/name.surname/nemo_4.2.0/cfgs/EQ_WAVES/EXP00'
# --------------------------------------------------------------------------- #
#orig_path='/home/PERSONALE/name.surname/nemo_4.2.0/cfgs/EQ_WAVES/EXP00'
orig_path=`pwd`
cd $orig_path
# create a tmp dir and move there
mkdir -p tmp_work
cd tmp_work

# --------------------------------------------------------------------------- #
# Copy the needed files to the target directory
# --------------------------------------------------------------------------- #
cp $orig_path/*.xml .
cp $orig_path/*namelist_* .
cp $orig_path/AGRIF* .
cp $orig_path/nemo .
cp $orig_path/namelist_cfg_0 namelist_cfg
# --------------------------------------------------------------------------- #
# Using the XIOS in server mode
# --------------------------------------------------------------------------- #
cp /home/software/terra/NEMO_env/xios/bin/xios_server.exe xios
# --------------------------------------------------------------------------- #
# Run the executable
# --------------------------------------------------------------------------- #
# xios attached 
# mpirun -np 4 ./nemo
# xios dettached 
mpirun -np 1 ./nemo : -np 1 xios
