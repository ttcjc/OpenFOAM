#!/bin/bash

#SBATCH --job-name=start
#SBATCH --output=start.out
#SBATCH --error=start.err
#SBATCH --time=6:00:00
#SBATCH --partition=compute
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3850
#SBATCH --account=su005-garmory21c


## Initialise OpenFOAM
module purge
module load GCC/11.3.0
module load flex/2.6.4
module load CMake/3.24.3
module load Qt5/5.15.5
module load OpenMPI/4.1.4
source ~/OpenFOAM/OpenFOAM-7/etc/bashrc


## Prepare Script
chmod +x Allstart

## Run
./Allstart
