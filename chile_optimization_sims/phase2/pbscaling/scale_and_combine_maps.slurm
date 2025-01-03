#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=16
#SBATCH --job-name=CMBS4_rescale_phase2
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=cpu
#SBATCH --account=mp107a

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=8

srun python -u scale_and_combine_maps.fullsky.py
# srun python -u scale_and_combine_maps.mflat.py
