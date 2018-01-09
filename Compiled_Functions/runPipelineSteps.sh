#!/bin/bash

step=0
for var in "$@"
do	
    if [[ $var == "cluster"* ]]; then
	cmd="-n 1 -P lightsheet ./run_compiled_matlab.sh ./run_distribute_parallelized.sh /usr/local/matlab-2017a/ $var sampleInput_$var.json 16 4"
    elif [[ $var == "local"* ]]; then
	cmd="-J \"$var\" -n 16 -P lightsheet ./run_compiled_matlab.sh ./run_localAP_fn.sh /usr/local/matlab-2017a/ sampleInput_$var.json"
    fi
    if [[ $step == 0 ]]; then
	cmd="bsub $cmd"
    else
	cmd="bsub -w \"done(${!step})\" -ti $cmd"
    fi
    $cmd
    ((step++))
done
