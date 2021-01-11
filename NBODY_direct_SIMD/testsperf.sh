#!/bin/bash

for nb_corps in '64' '128' '256' '512' '1024' '2048' '4096' '8192'
do
echo $nb_corps corps >> 'res406.txt'
./bin/NBODY_direct --in=../data/plummer_$nb_corps.nemo --tend=10.0 --dt=0.01 --soft=0.01 --sum >> 'res406.txt'
done


