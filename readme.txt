# requires a link to 'msms' in your bin !!!
# and the python library scipy + numpy installed

# Trying to estimate parameters of a selective sweep using ABC:
# Simply call the following line to simulate:
# 50 haploid individuals
# 1000 SNPs in the surveyed region of length rho = 4.N.r.L = 100 (N=effective size, r=per nucleotide recombination rate, L=length in bp)
# each combination of parameters are replicated 10 times
# the genomic region has a length of 100kb 
# msms used parameters read from the file: SAA_SaA_age_pos.prior
# Simulated models can be freely modified by adapting:
# line #24 in order to modify the msms command line
# the content of the prior file

# nInd=$1 # number of simulated individuals
# nSNP=$2 # number of SNPs to simulate
# rho=$3 # 4.N.r.L (N: number of individuals; r: recombination rate per bp; L: size of the genomic region)
# nRep=$4 # number of times a combination of parameters is replicated
# length=$5 # number of nucleotides of the surveyed region
# priorfile=$6 # file containing the combination of parameters. Each line of the priorfile will be replicated nRep times.

./ABCsweep.sh 50 1000 100 10 100000 SAA_SaA_age_pos.prior

# msms is used to simulate data
# msmscalc.py compute summary statistics from the msms output file, over the nRep replicates, and over the simulated region in a sliding window analysis
# it produces 3 output files:
# outputABC_positions.txt: list of positions of the slided windows
# outputABC_prior.txt: the used parameters. This file is redundant with the priorfile ... but have column names !
# outputABC_sumStats.txt:  the computed statistics per bin (average pi, std for pi within a bin, Tajima's D, Achaz's Y*, Pearson_r(position, pi) and its p_value

