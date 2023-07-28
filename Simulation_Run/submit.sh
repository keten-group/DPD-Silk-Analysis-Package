#!/bin/bash
#SBATCH --job-name="vf0.05"
#SBATCH -A p31412
#SBATCH -p long    ## partition
#SBATCH -N 2  ## number of nodes
#SBATCH -n 104  ## number of cores
#SBATCH --output=R-%x.%j.out
#SBATCH --ntasks-per-node=52  ## number of cores
#SBATCH -t 168:00:00

## cd $PBS_O_WORKDIR

# load necessary programs
module purge all
module load mpi

# run lammps
lmp=/projects/p31412/Diblock_DPD_Protocol_JWong_2017/lammps-17Nov16/src/lmp_mpi
mpirun -np 104 $lmp -in equil_shear_stretch.in -log lammps.log

## Provide Overwrite protection by sorting output files
# dir=./lammps_out
# if [ ! -d "$dir" ];then
# 	mpirun -np 28 $lmp -in equil_shear_stretch.in -log lammps.log
# 	mkdir $dir
# 	mv lammps.log *.dcd all_stress_11111.txt ./$dir
# else
# 	echo "Directory ${dir} already exists."
# fi