#!/bin/bash

#SBATCH --job-name=reconstruct
#SBATCH --output=reconstruct.out
#SBATCH --error=reconstruct.err
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64
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


## Load Parallel
module load GCCcore/10.3.0
module load parallel/20210622

## Prepare Script
chmod +x Allreconstruct

## Run
./Allreconstruct
