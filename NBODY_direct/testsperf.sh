#!/bin/bash

for nb_corps in '64' '128' '256' '512' '1024' '2048' '4096' '8192' '16384' '32768'
do
	echo $nb_corps corps >> 'SeqMutResTests.txt'
	./bin/NBODY_direct --in=../data/baredisk_$nb_corps-plummer_$nb_corps.nemo --tend=10.0 --dt=1.0 --soft=0.01 >> 'SeqMutResTests.txt'
done

