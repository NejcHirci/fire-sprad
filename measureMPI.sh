#!/bin/bash
#SBATCH --time=01:00:00      # walltime
#SBATCJ --reservation=fri    # reservation
#SBATCH --nodes=8            # number of nodes
#SBATCH --ntasks=8           # number of tasks
#SBATCH --cpus-per-task=4    # cpu-cores per task 
#SBATCH --hint=nomultithread # 1 thread per physical core 
#SBATCH --mem-per-cpu=500M   # memory per CPU core
#SBATCH -J "hyb_8_4"         # job name
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
index=1
H=500
W=500
EDGE_LEN=1
BLOCKS=1

while [ $index -le 5 ]
do
    myTimes=""
    printf "%d x %d | EDGE_LEN %d | BLOCKS %d \n" $H $W $EDGE_LEN $BLOCKS
    counter=1
    sumTime=0
    while [ $counter -le 10 ]
    do
        output=$(mpirun --map-by socket ./out $H $W $EDGE_LEN $BLOCKS);
        myTimes+="$output "
        printf "%f " $output
        sumTime=$(awk '{print $1+$2}' <<<"${sumTime} ${output}");
        counter=$(( counter + 1 ));
    done
    echo ""
    echo "---"
    mean=$(awk '{print $1/$2}' <<<"${sumTime} 10")
    echo "Mean"
    echo $mean
    echo "Standard deviation"
    echo $myTimes | awk -vM=$mean '{for(i=1;i<=NF;i++){sum+=($i-M)*($i-M)};print sqrt(sum/NF)}'
    echo "======================================="
    H=$(( H*2 ))
    W=$(( W*2 ))
    index=$(( index + 1))
done