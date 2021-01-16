#!/bin/bash
#SBATCH --time=48:00:00   # walltime
#SBATCH --ntasks=16  # number of processor cores (i.e. tasks)
#SBATCH --cpus-per-task=1 # number of cores per tak
#SBATCH --mem-per-cpu=3000M   # memory per CPU core
#SBATCH -J "pozar"   # job name
#SBATCH --output=result_mpi.txt

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun --map-by socket --mca btl self,vader,tcp ./out 
exit 1