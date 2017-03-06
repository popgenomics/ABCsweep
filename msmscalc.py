#!/usr/bin/python
# test: msms 10 10 -s tbs -SAA 200 -SaA tbs -SF 1e-2 -Smu 0
# cat prior.txt | msms 50 10 -s tbs -r 1 10000 -SAA tbs -SaA 100 -SF 0.01 -N 100000 >output_test.msms
import sys
inputFileName = sys.argv[1]
nIndiv = int(sys.argv[2])

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


