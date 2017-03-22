#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#define WIDTH 0.1
#define STEP 0.05

void window(float width, float step, std::vector <float> & minBin, std::vector <float> & maxBin);
void fillBins(std::vector<std::string> & bin_haplotypes, std::vector<std::string> const & haplotypes, std::vector<float> const & positions, std::vector<float> const & minBin, std::vector<float> const & maxBin, size_t const binID, const unsigned nIndiv);
void printBins(std::vector<std::string> const & bin_haplotypes);

void computePi(std::vector<std::string> const & bin_haplotypes, std::vector< std::vector< float> > & pi_bins, size_t const nIndiv, size_t const binID, const unsigned nCombIndiv);

float mean_over_replicates(std::vector<float> const & liste, std::vector<unsigned> const & valide_dataset, bool const neglectUnvalidated);


// arg1 = name of the msms outputfile
int main(int argc, char* argv[]){
	unsigned i(0);
	const std::string msmsFile(argv[1]); // mmsFile contains the name of msms output file to read
	const unsigned nIndiv = atoi(argv[2]); // number of individuals in the alignment
	const unsigned nCombParam = atoi(argv[3]); // number of combination of parameters
	const size_t nReplicate = atoi(argv[4]); // number of replicates per combination of parameters: n sim tot = nCombParam * nReplicate

	// define the boundaries of the surveyed windows over a sequence
	std::vector <float> minBin;
	std::vector <float> maxBin;
	window(WIDTH, STEP, minBin, maxBin);
	const size_t nBins(minBin.size());

	// number of combinations between nIndiv individuals
	const unsigned nCombIndiv(nIndiv * (nIndiv - 1) / 2);

	// open the msms output file
	std::ifstream fifo(msmsFile.c_str());

	if(fifo){ // if the file msmsFile exists
		std::string line;
		std::string word;
		unsigned int test(0);
		unsigned int test_bin(0);
		int nDataset(0); // count the number of simulated datasets over the whole msmsFile
		int nSNPs(0);
		unsigned int replicateID(0); // nDataset%nReplicate= {0, 1, ..., nReplicate-1, 0, 1, ..., nReplicate-1}
		size_t i = 0;
		size_t j = 0;
		size_t k = 0;
		int cntIndiv = -1;
		int cntSetOfParam = -1;

		std::vector <float> tmp_vector_float;

/*		unsigned valide_dataset[nReplicate];
		unsigned nSNPtotal[nReplicate];
		float pi_bins[nBins];
		std::string haplotypes[nIndiv];*/
		std::vector <unsigned> valide_dataset; // size = nReplicate
		std::vector <unsigned> nSNPtotal; // size = nReplicate
		std::vector < std::vector <float> > pi_bins; // size = [nBins][nReplicate]
		std::vector <std::string> haplotypes; // size = nIndiv haplotypes with all SNPs
		std::vector <std::string> bin_haplotypes; // size = nIndiv haplotypes with SNPs within bin
		std::vector <float> positions; // size = nSNPs
		while(std::getline(fifo, line)){ // read the msmsFile

			// 'segsites' line
			if(line[0]=='s'){ // if the line starts by a 's', then we expect: 'segsites: int'
				test=1;
				++nDataset;
				replicateID = (nDataset-1)%nReplicate; // the ID of the dataset over nReplicate replicates
				
				// if first dataset of the replicate
				if(replicateID == 0){
					++cntSetOfParam;
					valide_dataset.clear();
					nSNPtotal.clear();
					pi_bins.clear();

					// prepare pi_bins: pi_bins[bin][replicate]
					for(i=0; i<nBins; i++){
						pi_bins.push_back(tmp_vector_float);
					}

					std::cout << "set of param: " << cntSetOfParam << std::endl;
				}

				pi_bins.clear();
				haplotypes.clear();
				positions.clear();

				std::istringstream iss(line);
				i = 0;

				// read over the line containing 'segsites' to get the number of SNPs
				while( std::getline( iss, word, ' ') ){
					if( i==1 ){ // records the 2nd item of the line as the number of SNPs
						nSNPtotal.push_back(std::stoul(word));
					}
					++i;
				}


				// if no SNPs in the simulation
				if(nSNPtotal[replicateID] == 0){
					valide_dataset.push_back(0);
					// BALANCER ICI LES MISSING POUR LES VALEURS DE STATS DE CE REPLICAT
				}else{
					valide_dataset.push_back(1);
				}

				std::cout << "\treplicate_" << replicateID << "\tnSNP " << nSNPtotal[replicateID] << "\tvalidate" << valide_dataset[replicateID] << std::endl;


				continue; // stop reading the 'segites' line
			} // end of: "if the line contains the string: 'segsites'

		
			// 'positions' line	
			if(line[0]=='p'){ // if the line starts by a 'p', then we expect: 'positions: ...'
				cntIndiv = -1;
				i = 0; // count the number of elements in the line
				
				std::istringstream iss(line);

				// read over the line containing 'positions' to get the positions
				while( std::getline( iss, word, ' ') ){
					if( i>0 ){ // if word is not the first item of the line:
						positions.push_back(stof(word));
					}
					++i;
				}
		
				continue; // stop reading the 'positions' line
			} // end of: "if the line contains the string: 'positions'
			
			/* If the line is not empty, and doesn't start by '/', nor 'segsites', nor 'positions'*/
			//if(test == 1 && ligne!="" && ligne[0]!='/' && ligne[0]!='s' && ligne[0]!='p'){ 
			if(line!="" && line[0]!='/' && line[0]!='s' && line[0]!='p' && test==1){

				++cntIndiv;
				if(cntIndiv < nIndiv){
					haplotypes.push_back(line); // add sequences to <haplotypes>
				}

				// when all of the nIndiv were recorded within <haplotypes>
				if(cntIndiv == nIndiv-1){ // start of block of treatment of nIndiv haplotypes
					
					// loop over bins
					for(i=0; i<nBins; ++i){ // start loop over bins
						bin_haplotypes.clear(); // contains alignment of haplotypes for bin 'i'
						fillBins(bin_haplotypes, haplotypes, positions, minBin, maxBin, i, nIndiv);
				//		printBins(bin_haplotypes);

						computePi(bin_haplotypes, pi_bins, nIndiv, i, nCombIndiv);

					} // end of loop over bins


					// if the last replicate among replicates
					if(replicateID == (nReplicate - 1)){ // block for treatment of nReplicate
					// average over bins and replicates

					// clear
					for(i=0; i<nBins; ++i){
						pi_bins[i].clear();
					}
					
					} // end of block for treatment of nReplicate
				} // end of block of treatment of nIndiv haplotypes
			continue;
			}
			
		} // end of loop over the msmsFile
	}else{
		std::cerr << "ERROR: cannot oppen the file " << msmsFile << std::endl;
		exit(0);
	}	
	
	minBin.clear();
	maxBin.clear();
	return(0);	
}



void window(float width, float step, std::vector <float> & minBin, std::vector <float> & maxBin){
	float minB = 0.0;
	float maxB = width;

	minBin.push_back(minB);
	maxBin.push_back(maxB);

	while((minB + width) <= 1 && maxB < 1){
		minB += step;
		maxB += step;
		if(maxB > 1){
			maxB = 1;
		}
	minBin.push_back(minB);
	maxBin.push_back(maxB);
	}
}


void fillBins(std::vector<std::string> & bin_haplotypes, std::vector<std::string> const & haplotypes, std::vector<float> const & positions, std::vector<float> const & minBin, std::vector<float> const & maxBin, size_t const binID, const unsigned nIndiv){
	// function filling bins
//	std::cout << "scanning for SNPs in ]" << minBin[binID] << "-" << maxBin[binID] << "] interval" << std::endl;
	unsigned int test_bin = 0;
	size_t j(0);

	for(j=0; j<nIndiv; j++){
		bin_haplotypes.push_back("");
	}

	for(size_t i(0); i<positions.size(); ++i){
		if(test_bin == 1 && positions[i] > maxBin[binID]){
//			std::cout << "stop A: " << positions[i] << " greater than " << maxBin[binID] << std::endl;
			break;
		}
		if(positions[i] > minBin[binID] && positions[i] <= maxBin[binID]){
			test_bin = 1;
			for(j=0; j<nIndiv; j++){
				bin_haplotypes[j] += haplotypes[j][i];
			}
//			std::cout << positions[i] << " is in ]" << minBin[binID] << "-" << maxBin[binID] << "] interval" << std::endl;
			
		}
	}
//	std::cout << std::endl;
}


void printBins(std::vector<std::string> const & bin_haplotypes){
	size_t i{ 0 };
	std::cout << bin_haplotypes[0].size() << " SNPs" << std::endl;
	for(auto const & sequence : bin_haplotypes){
		std::cout << "sequence_" << ++i << " " << sequence << std::endl;
	}
}


void computePi(std::vector<std::string> const & bin_haplotypes, std::vector< std::vector< float> > & pi_bins, size_t const nIndiv, size_t const binID, const unsigned nCombIndiv){
	// compute Pi for the bin 'binID', recorded in pi_bins[binID][ vector of pi over replicates ]
	size_t const sizeBin(bin_haplotypes[0].size());
	float tmpPi(0.0);
	float meanPi(0.0);
	std::string allele="";

	// if no SNP in the bin --> return(pi = 0.0)
	if(sizeBin == 0){
		pi_bins[binID].push_back(0.0);
	}else{
	
	// if bins have SNPs:
		// array of pi of length equal to the number of SNPs
		float pi_over_SNPs[sizeBin];

		for(size_t i(0); i<sizeBin; ++i){ // loop over SNPs
			for(size_t j(0); j<nIndiv; ++j){ // loop over individuals
				allele = bin_haplotypes[j][i];
				pi_over_SNPs[i] += std::stoi(allele);
			} // end of loop over individuals

			tmpPi = (pi_over_SNPs[i] * (nIndiv - pi_over_SNPs[i])) / nCombIndiv; // tmpPi = #1 * #0 / nCombIndiv
			pi_over_SNPs[i] = tmpPi; // array of pi over SNPs --> use in the futur to compute pi_std and thetaW
			meanPi += tmpPi; // sum of pi values over SNPs
		} // end of loop over SNPs
		meanPi /= sizeBin; // meanPi = (sum of pi values / nCombIndiv) / number_of_SNP
		pi_bins[binID].push_back(meanPi);
	}
}


float mean_over_replicates(std::vector<float> const & liste, std::vector<unsigned> const & valide_dataset, bool const neglectUnvalidated){
	// return mean value of a statistics over replicates by rejecting unvalidated datasets (i.e: without SNP)
	float res(0.0);
	size_t cnt(0);

	for(size_t i(0); i<valide_dataset.size(); ++i){
		if(neglectUnvalidated == true){
			if(valide_dataset[i] == 1){
				++cnt;
				res += liste[i];
			}
		}else{
			++cnt;
			res+= liste[i];
		}
	}
	return(res/cnt);
}



