#!/bin/bash -l
# Use bash and pickup a basic login environment.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=defq             
#SBATCH --job-name=agents       
#SBATCH --output=array_%A-%a.log                
#SBATCH --array=1-40

pwd; hostname; date

echo "Running a program on $SLURM_JOB_NODELIST"

module load Workspace/v1

export JULIA_NUM_THREADS=4
~/julia-1.11.1/bin/julia agents.jl $SLURM_ARRAY_TASK_ID 3 3 --threads = 4
