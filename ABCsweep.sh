#!/bin/sh
#./ABCsweep.sh nInd nRep regionSize priorfile
nInd=$1 # number of simulated individuals
nSNP=$2 # number of SNPs to simulate
rho=$3 # 4.N.r.L (N: number of individuals; r: recombination rate per bp; L: size of the genomic region)
nRep=$4 # number of times a combination of parameters is replicated
length=$5 # number of nucleotides of the surveyed region 
priorfile=$6 # file containing the combination of parameters

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
cat ${priorfile}_tmp | msms $nInd $nSimul -s $nSNP -r $rho $length -SAA tbs -SaA tbs -SF tbs -Sp tbs -N 100000  >fifo &
./msmscalc.py fifo $nInd $nRep $nCombParam $length 
rm fifo
rm ${priorfile}_tmp
 
