#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -q iw-shared-6
#PBS -N surfaceVariation125
#PBS -o stdout
#PBS -e stderr 

cd $PBS_O_WORKDIR

module purge
#module load intel/15.0
#module load gcc/7.3.0
module load open64/4.5.1
module load mkl/11.2
module load mvapich2/2.1
module load git
module load anaconda3/2019.03
conda activate fenicsproject

python tap_sim.py
