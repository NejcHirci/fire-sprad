#!/bin/bash
#SBATCH --time=00:30:00      # walltime
#SBATCJ --reservation=fri    # reservation
#SBATCH --nodes=1            # number of nodes
#SBATCH --ntasks=4           # number of tasks
#SBATCH --cpus-per-task=1    # cpu-cores per task 
#SBATCH --hint=nomultithread # 1 thread per physical core 
#SBATCH --mem-per-cpu=500M   # memory per CPU core
#SBATCH -J "test"         # job name
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
H=2000
W=2000
EDGE_LEN=1
BLOCKS=0

mpirun --map-by socket ./out $H $W $EDGE_LEN $BLOCKS
exit 1