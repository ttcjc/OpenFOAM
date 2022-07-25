#!/bin/bash

#SBATCH --job-name=reconstruct
#SBATCH --output=reconstruct.out
#SBATCH --error=reconstruct.err
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=42
#SBATCH --mem-per-cpu=3850
#SBATCH --account=su005-garmory21c


## Load OpenFOAM
module purge
module load GCC/8.3.0
module load OpenMPI/3.1.4
module load OpenFOAM/7-20200508
module unload OpenMPI/3.1.4
module load openmpi/4.1.1
source $FOAM_BASH

## Load Parallel
module load GCCcore/10.3.0
module load parallel/20210622

## Prepare Script
chmod +x Allreconstruct

## Run
./Allreconstruct
