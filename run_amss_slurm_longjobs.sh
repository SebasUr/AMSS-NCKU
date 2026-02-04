#!/usr/bin/env bash
#SBATCH --job-name=AMSS
#SBATCH --partition=longjobs
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=1-12:00:00
#SBATCH --output=ERROUTS/slurm-%j.out
#SBATCH --error=ERROUTS/slurm-%j.err

# set -euo pipefail

# --- Modules ---
module purge
ml miniconda3/25.5.1
ml cuda/12.5
ml gnu14/14.2.0
ml openmpi5/5.0.7
ml libfabric/1.18.0

# --- Conda env ---
conda activate AMSS10

# --- Output directory naming ---
# This is picked up by AMSS_NCKU_Input.py
export AMSS_OUTPUT_DIR="GW150914_${SLURM_JOB_ID}_128OMP_1000"
export OMP_NUM_THREADS=4
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OMP_DISPLAY_ENV=VERBOSE

# Optional: keep OpenMP from spawning threads (code is MPI-first)
#export OMP_NUM_THREADS=1

# Run the Python workflow (compiles + runs ABE with mpirun -np 64)
# Using srun for Python is optional; leaving it as plain python is fine.
python AMSS_NCKU_Program.py
