Trying to estimate parameters of a selective sweep using ABC:
nInd=$1 # number of simulated individuals
nSNP=$2 # number of SNPs to simulate
rho=$3 # 4.N.r.L (N: number of individuals; r: recombination rate per bp; L: size of the genomic region)
nRep=$4 # number of times a combination of parameters is replicated
length=$5 # number of nucleotides of the surveyed region
priorfile=$6 # file containing the combination of parameters

./ABCsweep.sh 30 1000 800 10 100000 prior_SAA_SaA_age_pos.txt
./ABCsweep.sh 30 1000 200 5 100000 prior_SAA_SaA_age_pos.txt
