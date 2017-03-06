#!/bin/sh
#./ABCsweep.sh nInd nRep
nInd=$1
nRep=$2
priorfile=$3

# prepare replicates
nSimul=0
while IFS= read -r line; do
	for i in $(seq 1 1 $nRep); do
		nSimul=$(($nSimul+1))
		echo $line >> ${priorfile}_tmp
	done
done <$priorfile

mkfifo fifo
./msmscalc.py fifo $nInd $nRep &
cat ${priorfile}_tmp | msms $nInd $nSimul -s tbs -r 80 10000 -SAA tbs -SaA 100 -SF 0 -Sp 0.5 -N 100000  >fifo && rm fifo
rm ${priorfile}_tmp
 
