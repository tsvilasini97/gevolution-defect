#!/bin/bash
#BATCH --mail-user=tsvilasini97@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=N512_e11
#SBATCH --output=logs/N512_gevo.out
#SBATCH --error=logs/N512_gevo.err
#SBATCH --priority=50
#SBATCH --partition=shared-cpu
#SBATCH --ntasks=256
#SBATCH --mem-per-cpu=4000
#SBATCH --time=5:30:00
#SBATCH --parsable
#SBATCH --cpus-per-task=1

export OMP_NUM_THREADS=256

mpirun -np 256 ./gevolution -n 16 -m 16 -s settings.ini
