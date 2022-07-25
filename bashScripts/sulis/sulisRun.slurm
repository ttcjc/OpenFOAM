#!/bin/bash

#SBATCH --job-name=run
#SBATCH --output=run.out
#SBATCH --error=run.err
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
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

## Prepare Script
chmod +x Allrun

## Run
./Allrun
