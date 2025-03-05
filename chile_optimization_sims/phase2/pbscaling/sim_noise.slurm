#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=02:00:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=16
#SBATCH --job-name=CMBS4_phase2
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=cpu
#SBATCH --account=mp107a

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=16

# srun python -u sim_noise.lat.py
#srun python -u sim_noise.sat.py
srun python -u scale_and_combine_maps.fullsky.py

#srun python -u sim_noise.no_pbscaling.py
#srun python -u scale_and_combine_maps.no_pbscaling.py
