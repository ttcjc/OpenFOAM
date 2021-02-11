#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --job-name=meshConvert
#SBATCH --output=meshConvert.out
#SBATCH --error=meshConvert.err
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

## Prepare Script
chmod +x AllmeshConvert

## Run
./AllmeshConvert
