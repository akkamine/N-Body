#!/bin/bash

for nb_cors in '4' '8' '16' '32' '64'
do
	echo $nb_cors cors >> 'resMPIMutual.txt'
	mpirun -n $nb_cors -hostfile ~/hostfile --map-by node ./bin/NBODY_direct --in=../data/baredisk_32768-plummer_32768.nemo --tend=10.0 --dt=1.0 --soft=0.01 >> 'resMPIMutual.txt'

done
