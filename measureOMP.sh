#!/bin/bash
#SBATCH --time=02:30:00      # walltime
#SBATCH --reservation=fri
#SBATCH --nodes=1            # number of nodes
#SBATCH --ntasks=1           # number of tasks
#SBATCH --cpus-per-task=1    # cpu-cores per task 
#SBATCH --hint=nomultithread # 1 thread per physical core 
#SBATCH --mem-per-cpu=500M   # memory per CPU core
#SBATCH -J "omp_1"         # job name
#SBATCH --output=result_%x.txt

gcc -O2 -fopenmp -lm parallel_omp.c -o parallel_omp

# Define variables for fire spread
index=1
H=500
W=500
nthreads=$SLURM_CPUS_PER_TASK

printf "NTHREADS:%d\n" $nthreads

while [ $index -le 5 ]
do
    myTimes=""
    printf "%d x %d\n" $H $W
    counter=1
    sumTime=0
    while [ $counter -le 5 ]
    do
        output=$(srun -n1 --reservation=fri --cpus-per-task=$nthreads ./parallel_omp $H $W);
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