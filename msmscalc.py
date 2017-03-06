#!/usr/bin/python
# test: msms 10 10 -s tbs -SAA 200 -SaA tbs -SF 1e-2 -Smu 0
# cat prior.txt | msms 50 10 -s tbs -r 1 10000 -SAA tbs -SaA 100 -SF 0.01 -N 100000 >output_test.msms
import sys
inputFileName = sys.argv[1]
nIndiv = int(sys.argv[2])
nRep = int(sys.argv[3])

def getParams(x):
	res = {}
	x = x[3:].split("\t")
	for i in range(len(x)):
		if x[i][0]=="-":
			param = x[i][1:]
			value = float(x[i+1])
			res[param] = value
	return(res)

#nIndiv = 50
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


# parse the ms output file
#inputFileName = "output_test.msms"
totalData = parse_msms(inputFileName, nIndiv)


# compute the diversity for all SNPs
for i in totalData:
	if totalData[i]["segsites"] != 0:
		totalData[i]["pi"] = comp_pi(totalData[i]["haplotypes"], totalData[i]["segsites"])
		del totalData[i]["haplotypes"]	


# compute the average pi over replicates
bins = {}
for i in range(10):
	bins[i] = {}
	bins[i]["min"] = i/10.0
	bins[i]["max"] = (i+1)/10.0


pos_avg, pos_std = [], []
pi_avg, pi_std = [], []


bin_tmp = {}
for i in bins:
	bin_tmp[i] = {}
	bin_tmp[i]["positions"] = []
	bin_tmp[i]["pi"] = []

for i in range(nRep):
	for j in range(len(totalData[i]["positions"])):
		pos_tmp = totalData[i]["positions"][j]
		pi_tmp = totalData[i]["pi"][j]
		bin_value = [ k for k in bins if pos_tmp >= bins[k]["min"] and pos_tmp < bins[k]["max"] ]
		bin_tmp[bin_value]["positions"].append(pos_tmp)
		bin_tmp[bin_value]["pi"].append(pi_tmp)


# generates outputs
# for a fixed number of SNPs over all replicates
# parameters
header = ""
for i in totalData[0]["params"].keys():
	header += "{0}\t".format(i)
header += "\n"

simulations = ""
for i in totalData:
	if totalData[i]["segsites"] != 0:
		for j in totalData[i]["params"]:
			simulations += "{0}\t".format(totalData[i]["params"][j])
		simulations += "\n"
output = header + simulations
outfile = open("output_parameters.txt", "w")
outfile.write(output)


# diversity
simulations = ""
for i in totalData:
	if totalData[i]["segsites"] != 0:
		for j in totalData[i]["pi"]:
			simulations += "{0}\t".format(j)
		simulations += "\n"
output = simulations
outfile = open("output_pi.txt", "w")
outfile.write(output)

# positions
simulations = ""
for i in totalData:
	if totalData[i]["segsites"] != 0:
		for j in totalData[i]["positions"]:
			simulations += "{0}\t".format(j)
		simulations += "\n"
output = simulations
outfile = open("output_positions.txt", "w")
outfile.write(output)


