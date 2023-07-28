#!/bin/bash
#SBATCH --job-name="test"
#SBATCH -A p31412
#SBATCH -p short    ## partition
#SBATCH -N 1  ## number of nodes
#SBATCH -n 52  ## number of cores
#SBATCH --output=R-%x.%j.out
#SBATCH --ntasks-per-node=52  ## number of cores
#SBATCH -t 04:00:00

# load necessary programs
module purge all
module load mpi

# run lammps
lmp=/projects/p31412/Diblock_DPD_Protocol_JWong_2017/lammps-17Nov16/src/lmp_mpi # path to lammps executable
mpirun -np 52 $lmp -in equil_shear_stretch.in -log lammps.log