#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --job-name=mesh
#SBATCH --output=mesh.out
#SBATCH --error=mesh.err
#SBATCH --partition=compute
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=28
#SBATCH --account=a01-Garmory2020a
#SBATCH --mail-type=ALL
#SBATCH --mail-user=c.j.crickmore@lboro.ac.uk
#SBATCH --requeue

## Load OpenFOAM
module load site-local
module load openfoam/7.0/intel-2018.02
source $FOAM_BASH

## Prepare Scripts
chmod +x Allmesh

## Run
./Allmesh
