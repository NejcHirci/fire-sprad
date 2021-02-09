#!/bin/bash
#SBATCH --time=01:00:00      # walltime
#SBATCH --reservation=fri    # reservation
#SBATCH --nodes=3            # number of nodes
#SBATCH --ntasks-per-node=1  # number of tasks per node
#SBATCH --cpus-per-task=8    # cpu-cores per task for OMP
#SBATCH --hint=nomultithread # 1 thread per physical core 
#SBATCH --mem-per-cpu=500M   # memory per CPU core
#SBATCH -J "hyb_3_8"         # job name
#SBATCH --output=result_%x.txt

echo "SLURM_JOB_NNODES="$SLURM_JOB_NUM_NODES
echo "SLURM_JOB_CPUS_PER_NODE="$SLURM_JOB_CPUS_PER_NODE
echo "SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK
echo "SLURM_NTASKS="$SLURM_NTASKS
echo "SLURM_TASKS_PER_NODE="$SLURM_TASKS_PER_NODE
echo "---------------------"

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
    while [ $counter -le 5 ]
    do
        output=$(mpirun -x OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK -x OMP_PLACES=cores --bind-to socket ./out $H $W $EDGE_LEN $BLOCKS);
        myTimes+="$output "
        printf "%f " $output
        sumTime=$(awk '{print $1+$2}' <<<"${sumTime} ${output}");
        counter=$(( counter + 1 ));
    done
    echo ""
    echo "---"
    mean=$(awk '{print $1/$2}' <<<"${sumTime} 5")
    echo "Mean"
    echo $mean
    echo "Standard deviation"
    echo $myTimes | awk -vM=$mean '{for(i=1;i<=NF;i++){sum+=($i-M)*($i-M)};print sqrt(sum/NF)}'
    echo "======================================="
    H=$(( H*2 ))
    W=$(( W*2 ))
    index=$(( index + 1))
done