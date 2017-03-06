#!/bin/sh
#./ABCsweep.sh nInd nRep
nInd=$1 # number of simulated individuals
nRep=$2 # number of times a combination of parameters is replicated
length=$3 # number of nucleotides of the surveyed region 
priorfile=$4 # file containing the combination of parameters

# prepare replicates
nSimul=0
while IFS= read -r line; do
	for i in $(seq 1 1 $nRep); do
		nSimul=$(($nSimul+1))
		echo $line >> ${priorfile}_tmp
	done
done <$priorfile

nCombParam=$(wc -l $priorfile | awk '{print $1}')

mkfifo fifo
./msmscalc.py fifo $nInd $nRep $nCombParam $length&
cat ${priorfile}_tmp | msms $nInd $nSimul -s tbs -r 80 $length -SAA tbs -SaA 100 -SF 0 -Sp 0.5 -N 100000  >fifo
rm fifo
rm ${priorfile}_tmp
 
