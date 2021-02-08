#!/bin/bash
#SBATCH --time=03:30:00      # walltime
#SBATCJ --reservation=fri    # reservation
#SBATCH --nodes=1            # number of nodes
#SBATCH --ntasks=1           # number of tasks
#SBATCH --cpus-per-task=1    # cpu-cores per task 
#SBATCH --hint=nomultithread # 1 thread per physical core 
#SBATCH --mem-per-cpu=500M   # memory per CPU core
#SBATCH -J "serial"         # job name
#SBATCH --output=result_%x.txt

gcc -O2 -fopenmp -lm serial.c -o serial

index=1
H=500
W=500

while [ $index -le 5 ]
do
    myTimes=""
    printf "%d x %d\n" $H $W
    counter=1
    sumTime=0
    while [ $counter -le 10 ]
    do
        output=$(./serial $H $W);
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