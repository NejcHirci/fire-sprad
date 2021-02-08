#!/bin/bash
#SBATCH --time=00:10:00      # walltime
#SBATCH --reservation=fri    # fri
#SBATCH --nodes=4            # number of nodes
#SBATCH --ntasks=4           # number of tasks
#SBATCH --cpus-per-task=4    # cpu-cores per task 
#SBATCH --hint=nomultithread # 1 thread per physical core 
#SBATCH --mem-per-cpu=500M   # memory per CPU core
#SBATCH -J "test_hybrid_8_4"         # job name
#SBATCH --output=result_%x.txt

echo "SLURM_JOB_NNODES="$SLURM_JOB_NUM_NODES
echo "SLURM_JOB_CPUS_PER_NODE="$SLURM_JOB_CPUS_PER_NODE
echo "SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK
echo "SLURM_NTASKS="$SLURM_NTASKS
echo "SLURM_TASKS_PER_NODE="$SLURM_TASKS_PER_NODE
echo "---------------------"

# number of OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# OpenMP binding
export OMP_PLACES=cores

# Define variables for fire spread
H=8000
W=8000
EDGE_LEN=1
BLOCKS=1

mpirun --map-by socket ./out $H $W $EDGE_LEN $BLOCKS
exit 1