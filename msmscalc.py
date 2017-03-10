#!/usr/bin/python

# test: msms 10 10 -s tbs -SAA 200 -SaA tbs -SF 1e-2 -Smu 0
# cat prior.txt | msms 50 10 -s tbs -r 1 10000 -SAA tbs -SaA 100 -SF 0.01 -N 100000 >output_test.msms
import sys
import time
from numpy import nan
from numpy import mean
from numpy import nanmean
from numpy import std 
from numpy import sum as somme
from numpy import sqrt
from scipy.stats import pearsonr

minNumbSNP = 5 # statistics (pi, tajD, pearson's, etc ...) are only computed if a number of SNPs>= minNumbSNP is present within the studied fragment

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
	compteur = nCombParam * nRep
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
			if nLocus%(compteur/100) == 0:
				print("simulation {0} over a total of {1}: {2}".format(nLocus, compteur, time.strftime("%H:%M:%S")))
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
#			print("{0}: {1}".format(nIndTmp, i))
			res[nLocus]["haplotypes"].append(i)
			if nIndTmp == (nIndiv-1):
#				print("{0}\n{1}".format(nLocus, "\n".join(res[nLocus]["haplotypes"])))
				tmp = comp_differences(res[nLocus]["haplotypes"], res[nLocus]["segsites"])
#				print(a["kxy_allSNP"])
				res[nLocus]['kxy'] = tmp['kxy_allSNP'] 
				res[nLocus]['kxy_singletons'] = tmp['kxy_singletons'] 
				del res[nLocus]['haplotypes']	
	infile.close()
	return(res)


def comp_differences(alignement, nSegSites):
	n = len(alignement) # number of individuals in the alignement
	res = {}
	res["kxy_allSNP"] = [] # vector containing pi for different SNPs, the probability to sample 2 different alleles
	res["kxy_singletons"] = [] # vector containing pi for different SNPs for singletons only
	for i in range(nSegSites):
		n0 = 0
		n1 = 0
		for j in alignement:
			if j[i] == '0':
				n0 += 1
			if j[i] == '1':
				n1 += 1
		pi = n0*n1
		res["kxy_allSNP"].append(pi)
		if n0==1 or n1==1: # if only one derived allele or one ancestral allele
			res["kxy_singletons"].append(pi)
	return(res)


def comp_pi(alignement, nSegSites):
	n = len(alignement) # number of individuals in the alignement
	res = {}
	res["pi_allSNP"] = [] # vector containing pi for different SNPs, the probability to sample 2 different alleles
	res["pi_singletons"] = [] # vector containing pi for different SNPs for singletons only
	for i in range(nSegSites):
		n0 = 0
		n1 = 0
		nComparaisons = n*(n-1.0) / 2.0
		for j in alignement:
			if j[i] == '0':
				n0 += 1
			if j[i] == '1':
				n1 += 1
		pi = n0*n1/nComparaisons
		res["pi_allSNP"].append(pi)
		if n0==1 or n1==1: # if only one derived allele or one ancestral allele
			res["pi_singletons"].append(pi)
	res["nSingletons"] = len(res["pi_singletons"])
	res["nSNP"] = len(res["pi_allSNP"])
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


def tajimaD(nInd, pi, nS):
	# nInd = number of individuals in the alignment
	# pi = sum(Kxy over SNPs and pairwise comparisons) / number_of_pairwise_comparisons
	# nS = number of SNPs within the alignment
	# a1 and a2
	minNumbSNP = 5
	if nS < minNumbSNP:
		return(nan)
	else:
		a1, a2 = 0.0, 0.0
		for i in range(nInd-1):
			a1 += 1.0/(i+1.0)
			a2 += 1.0/((i+1)**2)
		# b1 and b2
		b1 = (nInd + 1.0) / (3.0 * (nInd - 1.0))
		b2 = 2.0 * (nInd**2 + nInd + 3.0) / (9.0*nInd * (nInd - 1.0))
		# c1 and c2
		c1 = b1 - 1.0 / a1
		c2 = b2 - (nInd + 2.0) / (a1 * nInd) + a2/(a1**2.0)
		# e1 and e2
		e1 = c1/a1
		e2 = c2/(a1**2 + a2)
		# pi is assumed to be: sum(pi over SNPs)/size, let's compute thetaW
		thetaW = nS / a1
		# denominateur
		denominateur = e1*nS + e2*nS*(nS - 1.0)
		denominateur = sqrt(denominateur)
		# test
	#	print("a1 = {0}\na2 = {1}\nb1 = {2}\nb2 = {3}\nc1 = {4}\nc2 = {5}\ne1 = {6}\ne2 = {7}\npi = {8}\nthetaW = {9}\ndenoM = {10}".format(a1, a2, b1, b2, c1, c2, e1, e2, pi, thetaW, denominateur))
		#tajima D
		return((pi - thetaW) / denominateur)


def achazY(n, kxy):
	# n = number of individuals
	# kx = vector of pi values for each singletons over the bin
	minNumbSNP = 5
	nPairwiseComp = n * (n-1.0)/2.0
	singletons = [ i for i in kxy if i==(n-1) ]
	nSingletons = len(kxy)
	if nSingletons < minNumbSNP: # if less than minNumbSNP of SNPs are present in the alignement --> return "nan" value
		return(nan)
	else:
		an = 0.0
		bn = 0.0
		for i in range(n-1):
			an += 1.0/(i+1)
			bn += 1.0/(i+1)**2
		# f*
		f = (n - 3.0) / (an * (n-1) - n)
		# alpha
		alpha = f**2 * (an - n / (n-1.0)) + f * (an * (4.0*(n+1)/((n-1.0)**2)) - 2 * ((n+3.0)/(n-1))) - an*(8*(n+1.0)/(n*(n-1)**2)) + (n**2 + n + 60.0)/(3*n * (n-1))
		# beta
		beta = f**2 * (bn - (2*n - 1.0)/(n-1.0)**2) + f * (bn * 8.0/(n-1.0) - an * 4.0/(n*(n-1.0)) - (n**3 + 12*n**2 - 35*n + 18.0)/(n*(n-1.0)**2)) - bn * 16.0/(n*(n-1)) + an * 8.0 / (n**2 * (n-1.0)) + 2 * (n**4 + 110.0*n**2 - 255*n + 126.0) / (9.0 * n**2 * (n-1)**2)
		# mu
		mu = an - n/(n-1.0)
		# theta
		theta = (1.0 * nSingletons) / mu
		# pi singleton: sum of kxy for singletons only
		pi_singleton = 0.0
		for i in singletons:
			pi_singleton += i
		pi_singleton /= nPairwiseComp
		# compute Y* from Achaz
		numerateur = pi_singleton - f*nSingletons/an
		denominateur = alpha * theta + beta * theta**2
		return((1.0 * numerateur)/(1.0 * denominateur))


def calc_window(rep, nRep, bins, regionSize, nIndiv):
	# function used to return list of 'positions' and 'associated pi' within different bins, over replicates in the
	# rep = the id of the surveyed replicated simulations
	# nRep = number of times the simulation had been replicated
	nPairwiseComp = (nIndiv * (nIndiv - 1.0)) / 2.0
	bins_tmp = {} # bins_tmp[bin_Id][replicate_Id][['positions'], ['kxy']]
	for i in bins:
		bins_tmp[i] = {}
		for j in range(nRep):
			bins_tmp[i][j] = {} # bins_tmp[bin ID][rep ID]
			bins_tmp[i][j]['positions'] = [] # bins_tmp[bin ID][rep ID]
	#		bins_tmp[i][j]['pi'] = [] # bins_tmp[bin ID][rep ID]
			bins_tmp[i][j]['kxy'] = [] # bins_tmp[bin ID][rep ID]
	for i in range(nRep): # loop over the replicates
		for j in range(len(totalData[rep*nRep + i]['positions'])): # loop over positions
			pos_tmp = totalData[rep*nRep + i]['positions'][j]
			kxy_tmp = totalData[rep*nRep + i]['kxy'][j] # sum of pairwise differences for a SNP
			bin_value = [ k for k in bins if pos_tmp >= bins[k]['min'] and pos_tmp < bins[k]['max'] ]
			for l in bin_value:
				bins_tmp[l][i]['positions'].append(pos_tmp)
				bins_tmp[l][i]['kxy'].append(kxy_tmp)
	#		print(" ".join([ str(totalData[rep*nRep]['params'][k]) for k in totalData[rep*nRep]['params'] ])) # print test of homogeneity in parameters over replicates
	res = {}
	for i in bins: # compute statistics #1 over bins and over #2 replicates
		meanPi_tmp = []
		stdPi_tmp = []
		pearsonR_tmp = []
		pearsonPval_tmp = []
		tajimaD_tmp = []
		achazY_tmp = []
		pos_tmp = []
		size_tmp = regionSize * (bins[i]['max'] - bins[i]['min'])
		for j in range(nRep):
			nS = len(bins_tmp[i][j]['positions']) # number of SNPs
	#		meanPi_tmp.append(mean(bins_tmp[i][j]['pi']) / size_tmp)
			meanPi_tmp.append(sum(bins_tmp[i][j]['kxy']) /(1.0 * nPairwiseComp)) # sum(Kxy over comparisons at a SNP) / number_pairwise_comparison 
			stdPi_tmp.append(stdCustom(bins_tmp[i][j]['kxy'], meanPi_tmp[j], regionSize))
			if len(bins_tmp[i][j]['positions'])<minNumbSNP:
				pearson = [nan, nan]
			else:
				pearson = pearsonr(bins_tmp[i][j]['positions'], bins_tmp[i][j]['kxy'])
			pearsonR_tmp.append(pearson[0])
			pearsonPval_tmp.append(pearson[1])
			tajimaD_tmp.append(tajimaD(nIndiv, meanPi_tmp[j], nS))
			achazY_tmp.append(achazY(nIndiv, bins_tmp[i][j]['kxy']))
			pos_tmp.append(mean(bins_tmp[i][j]['positions']))
			#def stdCustom(liste, moyenne, longueurRegion):
		res[i] = {}
		res[i]['pi_avg'] = round(nanmean(meanPi_tmp)/size_tmp, 5)
		res[i]['pi_std'] = round(nanmean(stdPi_tmp)/size_tmp, 5)
		res[i]['pearson_r'] = round(nanmean(pearsonR_tmp), 5)
		res[i]['pearson_pval'] = round(nanmean(pearsonPval_tmp), 5)
		res[i]['tajD'] = round(nanmean(tajimaD_tmp), 5)
		res[i]['achazY'] = round(nanmean(achazY_tmp), 5)
		res[i]['position'] = round(nanmean(pos_tmp), 5)
	return(res)


# parse the ms output file
#inputFileName = "output_test.msms"
totalData = parse_msms(inputFileName, nIndiv)


# compute the diversity for all SNPs, for all alignments
# commented block because: loop is now imbedded in the parse_msms function to save memory
# now: the alignments are treated during reading, not only after
#for i in totalData:
#	if totalData[i]['segsites'] != 0:
#		#tmp = comp_pi(totalData[i]['haplotypes'], totalData[i]['segsites']) # first version: return the average pi
#		#totalData[i]['pi'] = tmp['pi_allSNP'] 
#		#totalData[i]['pi_singletons'] = tmp['pi_singletons'] 
#		tmp = comp_differences(totalData[i]['haplotypes'], totalData[i]['segsites']) # second version: return the sum of differences
#		totalData[i]['kxy'] = tmp['kxy_allSNP'] 
#		totalData[i]['kxy_singletons'] = tmp['kxy_singletons'] 
#		del totalData[i]['haplotypes']	


# compute the average pi over replicates
width = 0.1
step = 0.05

# define the bin boundaries
bins = window(width, step)


# prepare the output files:
stats_out = ""
for i in bins:
	stats_out += 'pi_avg_bin{0}\t'.format(i)
	stats_out += 'pi_std_bin{0}\t'.format(i)
	stats_out += 'tajD_bin{0}\t'.format(i)
	stats_out += 'achazY_bin{0}\t'.format(i)
	stats_out += 'pearson_r_bin{0}\t'.format(i)
	stats_out += 'pearson_pval_bin{0}\t'.format(i)
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
	print("simulation {0} over a total of {1}: {2}".format(i, nCombParam, time.strftime("%H:%M:%S")))
	# TODO: calling this function tajimaD(nInd, pi, nS, size)
	a = calc_window(i, nRep=nRep, bins=bins, regionSize=regionSize, nIndiv=nIndiv) # get summary statistics per bin for the replicate i
	for j in bins:
		stats_out += "{0}\t".format(a[j]['pi_avg'])
		stats_out += "{0}\t".format(a[j]['pi_std'])
		stats_out += "{0}\t".format(a[j]['tajD'])
		stats_out += "{0}\t".format(a[j]['achazY'])
		stats_out += "{0}\t".format(a[j]['pearson_r'])
		stats_out += "{0}\t".format(a[j]['pearson_pval'])
		pos_out += "{0}\t".format(a[j]['position'])
	stats_out = stats_out.strip() + "\n"
	pos_out = pos_out.strip() + "\n"


outfile = open("outputABC_prior.txt", "w")
outfile.write(params_out)
outfile.close()

outfile = open("outputABC_positions.txt", "w")
outfile.write(pos_out)
outfile.close()

outfile = open("outputABC_sumStats.txt", "w")
outfile.write(stats_out)
outfile.close()


