#!/bin/sh
#./ABCsweep.sh nInd nRep regionSize priorfile
nInd=$1 # number of simulated individuals
nSNP=$2 # number of SNPs to simulate
rho=$3 # 4.N.r.L (N: number of individuals; r: recombination rate per bp; L: size of the genomic region)
length=$4 # number of nucleotides of the surveyed region 


mkfifo fifo
msms $nInd 2500 -s $nSNP -r $rho $length -N 100000  >fifo &
./msmscalc.py fifo $nInd 1 2500 $length 
rm fifo
 
 
