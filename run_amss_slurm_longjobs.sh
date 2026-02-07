#!/usr/bin/env bash
#SBATCH --job-name=AMSS_OPT
#SBATCH --partition=longjobs
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=1-12:00:00
#SBATCH --output=ERROUTS/slurm-%j.out
#SBATCH --error=ERROUTS/slurm-%j.err

module purge
ml miniconda3/25.5.1
ml cuda/12.5
ml gnu14/14.2.0
ml openmpi5/5.0.7
ml libfabric/1.18.0

conda activate AMSS10

export AMSS_OUTPUT_DIR="GW150914_${SLURM_JOB_ID}_TEST"
export OMP_NUM_THREADS=1

python AMSS_NCKU_Program.py
