#!/usr/bin/python

# test: msms 10 10 -s tbs -SAA 200 -SaA tbs -SF 1e-2 -Smu 0
# cat prior.txt | msms 50 10 -s tbs -r 1 10000 -SAA tbs -SaA 100 -SF 0.01 -N 100000 >output_test.msms
import sys
from numpy import mean
from numpy import std 
from numpy import sum as somme
from numpy import sqrt
from scipy.stats import pearsonr


inputFileName = sys.argv[1]
nIndiv = int(sys.argv[2])
nRep = int(sys.argv[3])
nCombParam = int(sys.argv[4])
regionSize = int(sys.argv[5])


def stdCustom(liste, moyenne, longueurRegion):
	# 'moyenne' = mean of ['liste' (size L) + longueurRegion x 0]
	res = 0.0
	for i in liste:
		res += (i - moyenne)**2
	res /= longueurRegion
	return(sqrt(res))


def getParams(x):
	res = {}
	x = x[3:].split("\t")
	for i in range(len(x)):
		if x[i][0]=="-":
			param = x[i][1:]
			value = float(x[i+1])
			res[param] = value
	return(res)


def parse_msms(x, nIndiv):
	# returns the simulated alignments
	# x = msms output file
	res = {}
	infile = open(x, "r")
	nLocus = -1
	nIndTmp = -1
	test = -9999
	for i in infile:
		i = i.strip()
		if "//" in i:
			test = 1
			nLocus += 1
			res[nLocus] = {}
			if i != "//":
				res[nLocus]["params"] = getParams(i)
			continue
		if "segsites" in i and nLocus >= 0:
			res[nLocus]["segsites"] = int(i.split(" ")[1])
			if res[nLocus]["segsites"] == 0:
				continue
			continue
		if "positions" in i and nLocus >= 0:
			res[nLocus]["positions"] = [ float(pos) for pos in i.split(" ")[1::] ]
			res[nLocus]["haplotypes"] = []
			nIndTmp = -1
			continue
		if i != "" and test == 1:
			nIndTmp += 1
			#print("{0}: {1}".format(nIndTmp, i))
			res[nLocus]["haplotypes"].append(i)
	infile.close()
	return(res)


def comp_pi(alignement, nSegSites):
	res = []
	for i in range(nSegSites):
		n0 = 0
		n1 = 0
		nComparaisons = len(alignement)*(len(alignement)-1) / 2.0
		for j in alignement:
			if j[i] == '0':
				n0 += 1
			if j[i] == '1':
				n1 += 1
		res.append(n0*n1/nComparaisons)
	return(res)


def window(width, step):
	bins = {}
	cnt = 0
	minB = 0
	maxB = width
	bins[cnt] = {}
	bins[cnt]["min"] = minB 
	bins[cnt]["max"] = maxB 
	while(round((minB+width), 5) <= 1 and maxB < 1):
		cnt += 1
		minB += step
		maxB += step
		if maxB > 1:
			maxB = 1
		bins[cnt] = {}
		bins[cnt]["min"] = round(minB, 5)
		bins[cnt]["max"] = round(maxB, 5)
	return(bins)


def calc_window(rep, nRep, bins):
	# function used to return list of 'positions' and 'associated pi' within different bins, over replicates in the
	# rep = the id of the surveyed replicated simulations
	# nRep = number of times the simulation had been replicated
	bins_tmp = {}
	for i in bins:
		bins_tmp[i] = {}
		bins_tmp[i]["positions"] = []
		bins_tmp[i]["pi"] = []
	for i in range(nRep): # loop over the replicates
		for j in range(len(totalData[rep*nRep + i]["positions"])): # loop over positions
			pos_tmp = totalData[rep*nRep + i]["positions"][j]
			pi_tmp = totalData[rep*nRep + i]["pi"][j]
			bin_value = [ k for k in bins if pos_tmp >= bins[k]["min"] and pos_tmp < bins[k]["max"] ]
			for l in bin_value:
				bins_tmp[l]["positions"].append(pos_tmp)
				bins_tmp[l]["pi"].append(pi_tmp)
#		print(" ".join([ str(totalData[rep*nRep]['params'][k]) for k in totalData[rep*nRep]['params'] ])) # print test of homogeneity in parameters over replicates
	return(bins_tmp)


# parse the ms output file
#inputFileName = "output_test.msms"
totalData = parse_msms(inputFileName, nIndiv)


# compute the diversity for all SNPs
for i in totalData:
	if totalData[i]["segsites"] != 0:
		totalData[i]["pi"] = comp_pi(totalData[i]["haplotypes"], totalData[i]["segsites"])
		del totalData[i]["haplotypes"]	


# compute the average pi over replicates
width = 0.1
step = 0.05

bins = window(width, step)


# prepare the outputs:
stats_out = ""
for i in bins:
	stats_out += "pi_avg_bin_{0}\t".format(i)


for i in bins:
	stats_out += "pi_std_bin_{0}\t".format(i)


for i in bins:
	stats_out += "pearson_r_bin_{0}\tpearson_pval_bin_{0}\t".format(i)
stats_out = stats_out.strip() + "\n"


pos_out = ""
for i in bins:
	pos_out += "SNP_{0}\t".format(i)
pos_out = pos_out.strip() + "\n"


params_out = ""
for i in totalData[0]["params"].keys():
	params_out += "{0}\t".format(i)
params_out += "\n"


for i in range(nCombParam):
	for j in totalData[i*nRep]["params"]:
		 params_out += "{0}\t".format(totalData[i*nRep]["params"][j])
	params_out = params_out.strip() + "\n"


# treat the replicated datasets
for i in range(nCombParam): # loop over combination of parameters
	a = calc_window(i, nRep=nRep, bins=bins) # to pool bins of replicated simulations
	res_tmp = ""
	for j in bins:
		# pi_avg = sum(pairwise_differences) / (#replicates * regionSize)
		L_tmp = nRep*regionSize*(max(a[j]['positions'])-min(a[j]['positions'])) # L_tmp = #replicates x regionSize x relative_bin_size
		mean_tmp = somme(a[j]['pi'])/L_tmp # mean_pi = sum of pi over SNPs and over replicates / L_tmp
		res_tmp += "{0}\t".format(round(mean_tmp, 5))
		pos_out += "{0}\t".format(round(mean(a[j]['positions']), 5))
	pos_out = pos_out.strip() + "\n"
	for j in bins:
		std_tmp = stdCustom(a[j]['pi'], mean_tmp, L_tmp)
#		 res_tmp += "{0}\t".format(round(std(a[j]['pi']), 5))
		res_tmp += "{0}\t".format(round(std_tmp, 5))
	for j in bins:
		tmp = pearsonr(a[j]['positions'], a[j]['pi'])
		res_tmp += "{0}\t{1}\t".format(round(tmp[0], 5), round(tmp[1], 5))
	res_tmp = res_tmp.strip() + "\n"
	stats_out += res_tmp


outfile = open("outputABC_prior.txt", "w")
outfile.write(params_out)
outfile.close()

outfile = open("outputABC_positions.txt", "w")
outfile.write(pos_out)
outfile.close()

outfile = open("outputABC_sumStats.txt", "w")
outfile.write(stats_out)
outfile.close()



