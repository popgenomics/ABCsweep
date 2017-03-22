#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#define WIDTH 0.1
#define STEP 0.05

void window(float width, float step, std::vector <float> & minBin, std::vector <float> & maxBin);
void fillBins(std::vector<std::string> & bin_haplotypes, std::vector<std::string> const & haplotypes, std::vector<float> const & positions, std::vector<float> const & minBin, std::vector<float> const & maxBin, size_t const binID);

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


/*		unsigned valide_dataset[nReplicate];
		unsigned nSNPtotal[nReplicate];
		float pi_bins[nBins];
		std::string haplotypes[nIndiv];*/
		std::vector <unsigned> valide_dataset; // size = nReplicate
		std::vector <unsigned> nSNPtotal; // size = nReplicate
		std::vector <float> pi_bins; // size = nBins
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
					valide_dataset.clear();
					nSNPtotal.clear();
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
				if(valide_dataset[replicateID] == 1){
					++cntIndiv;
					if(cntIndiv < nIndiv){
						haplotypes.push_back(line); // add sequences to <haplotypes>
					}

					// when all of the nIndiv were recorded within <haplotypes>
					if(cntIndiv == nIndiv-1){ // start of block of treatment of nIndiv haplotypes
						
						// loop over bins
						for(i=0; i<nBins; ++i){ // start loop over bins
							bin_haplotypes.clear();
							fillBins(bin_haplotypes, haplotypes, positions, minBin, maxBin, i);

						} // end of loop over bins

	
						// if the last replicate among replicates
						if(replicateID == (nReplicate - 1)){ // block for treatment of nReplicate
							// do averages over replicates
							for(i=0; i<nIndiv; ++i){ // loop over replicates after recorded the last one
								
							}// end of loop over replicates after recorded the last one
						} // end of block for treatment of nReplicate
					} // end of block of treatment of nIndiv haplotypes
				}
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


void fillBins(std::vector<std::string> & bin_haplotypes, std::vector<std::string> const & haplotypes, std::vector<float> const & positions, std::vector<float> const & minBin, std::vector<float> const & maxBin, size_t const binID){
	std::cout << "scanning for SNPs in ]" << minBin[binID] << "-" << maxBin[binID] << "] interval" << std::endl;
	unsigned int test_bin = 0;
	for(size_t i(0); i<positions.size(); ++i){
		if(test_bin == 1 && positions[i] > maxBin[binID]){
			std::cout << "stop A: " << positions[i] << " greater than " << maxBin[binID] << std::endl;
			break;
		}
		if(positions[i] > minBin[binID] && positions[i] <= maxBin[binID]){
			test_bin = 1;
			std::cout << positions[i] << " is in ]" << minBin[binID] << "-" << maxBin[binID] << "] interval" << std::endl;
		}
	}
	std::cout << std::endl;
}

