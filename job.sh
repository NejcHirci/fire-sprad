#!/bin/bash
#SBATCH --time=00:30:00   # walltime
#SBATCH --nodes=2   # number of nodes
#SBATCH --ntasks=2  # number of processor cores (i.e. tasks)
#SBATCH --cpus-per-task=2 # number of nodes per task
#SBATCH --mem-per-cpu=3000M   # memory per CPU core
#SBATCH -J "pozar"   # job name
#SBATCH --output=result_%j.txt

echo "SLURM_JOB_NNODES="$SLURM_JOB_NUM_NODES
echo "SLURM_JOB_CPUS_PER_NODE="$SLURM_JOB_CPUS_PER_NODE
echo "SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK
echo "SLURM_NTASKS="$SLURM_NTASKS
echo "SLURM_TASKS_PER_NODE="$SLURM_TASKS_PER_NODE
echo "---------------------"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun --map-by core --mca btl self,vader,tcp ./out 
exit 1