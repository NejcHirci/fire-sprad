#!/bin/bash
#SBATCH --time=00:30:00      # walltime
#SBATCH --reservation=fri    # reservation
#SBATCH --nodes=3            # number of nodes
#SBATCH --ntasks-per-node=1  # always use one task per node
#SBATCH --cpus-per-task=8    # cpu-cores per task 
#SBATCH --hint=nomultithread # 1 thread per physical core 
#SBATCH --mem-per-cpu=1000M   # memory per CPU core
#SBATCH -J "test_hybd_3_8_test" # job name
#SBATCH --output=result_%x.txt

echo "SLURM_JOB_NNODES="$SLURM_JOB_NUM_NODES
echo "SLURM_JOB_CPUS_PER_NODE="$SLURM_JOB_CPUS_PER_NODE
echo "SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK
echo "SLURM_NTASKS="$SLURM_NTASKS
echo "SLURM_TASKS_PER_NODE="$SLURM_TASKS_PER_NODE
echo "---------------------"

# Define variables for fire spread
H=8000
W=8000
EDGE_LEN=1
BLOCKS=1

mpirun -x OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK -x OMP_PLACES=cores --bind-to socket ./out $H $W $EDGE_LEN $BLOCKS
exit 1