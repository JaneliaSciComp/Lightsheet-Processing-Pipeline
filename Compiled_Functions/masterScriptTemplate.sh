#!/bin/bash
bsub -n 1 -P lightsheet ./run_compiled_matlab.sh ./run_distribute_parallelized.sh /usr/local/matlab-2017a/ clusterPT sampleInput_clusterPT.json 16 4
bsub -w "done(clusterPT)" -ti -n 1 -P lightsheet ./run_compiled_matlab.sh ./run_distribute_parallelized.sh /usr/local/matlab-2017a/ clusterMF sampleInput_clusterMF.json 16 4
bsub -w "done(clusterMF)" -J "localAP" -ti -n 16 -P lightsheet ./run_compiled_matlab.sh ./run_localAP_fn.sh /usr/local/matlab-2017a/ sampleInput_localAP.json
bsub -w "done(localAP)" -ti -n 1 -P lightsheet ./run_compiled_matlab.sh ./run_distribute_parallelized.sh /usr/local/matlab-2017a/ clusterTF sampleInput_clusterTF.json 16 4
bsub -w "done(clusterTF)" -J "localEC" -ti -n 16 -P lightsheet ./run_compiled_matlab.sh ./run_localEC_fn.sh /usr/local/matlab-2017a/ sampleInput_localEC.json
bsub -w "done(localEC)" -ti -n 1 -P lightsheet ./run_compiled_matlab.sh ./run_distribute_parallelized.sh /usr/local/matlab-2017a/ clusterCS sampleInput_clusterCS.json 16 4
bsub -w "done(clusterCS)" -ti -n 1 -P lightsheet ./run_compiled_matlab.sh ./run_distribute_parallelized.sh /usr/local/matlab-2017a/ clusterFR sampleInput_clusterFR.json 16 4


