#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --job-name=reconstruct
#SBATCH --output=reconstruct.out
#SBATCH --error=reconstruct.err
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=a01-Garmory2020a
#SBATCH --mail-type=ALL
#SBATCH --mail-user=c.j.crickmore@lboro.ac.uk
#SBATCH --requeue

## Load OpenFOAM
module load site-local
module load openfoam/7.0/intel-2018.02
source $FOAM_BASH

## Prepare Scripts
chmod +x Allreconstruct

## Run
./Allreconstruct
