#!/bin/sh
#./ABCsweep.sh nInd nRep
nInd=$1
nRep=$2
mkfifo tmp_file
./msmscalc.py tmp_file $nInd &
cat prior_v2.txt | msms $nInd $nRep -s tbs -r 80 10000 -SAA tbs -SaA 100 -SF 0 -Sp 0.5 -N 100000  >tmp_file && rm tmp_file
 
