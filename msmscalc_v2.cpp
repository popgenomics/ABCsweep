#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#define WIDTH 0.1
#define STEP 0.05

// compil with: g++ msmscalc_v2.cpp -std=c++17 -O3 -o msmscalc

void window(float width, float step, std::vector <float> & minBin, std::vector <float> & maxBin);
void fillBins(std::vector<std::string> & bin_haplotypes, std::vector<std::string> const & haplotypes, std::vector<float> const & positions, std::vector<float> const & minBin, std::vector<float> const & maxBin, size_t const binID, const unsigned nIndiv);
void printBins(std::vector<std::string> const & bin_haplotypes);

void computePi(std::vector<std::string> const & bin_haplotypes, std::vector< std::vector< float> > & pi_avg_bins, std::vector< std::vector< float> > & pi_std_bins, std::vector< std::vector< float> > & thetaW_bins, std::vector< std::vector< float> > & tajimaD_bins, size_t const nIndiv, size_t const binID, const unsigned nCombIndiv, unsigned int replicateID, float const an, float const an2);

float tajimaD(unsigned const nIndiv, float const pi, float const thetaW, float const an, float an2, const unsigned nS);

float mean_piStats(std::vector<float> const & liste);

void write_header(std::string const outputFile, size_t const nBins);
void write_results(std::string const outputFile, float const pi_avg_bin[], float const pi_std_bin[], float const thetaW_bin[], float const tajimaD_bin[], size_t const nBins);

// arg1 = name of the msms outputfile
int main(int argc, char* argv[]){
	unsigned i(0);
	const std::string msmsFile(argv[1]); // mmsFile contains the name of msms output file to read
	const unsigned nIndiv = atoi(argv[2]); // number of individuals in the alignment
	const unsigned nCombParam = atoi(argv[3]); // number of combination of parameters
	const size_t nReplicate = atoi(argv[4]); // number of replicates per combination of parameters: n sim tot = nCombParam * nReplicate
	const std::string outputFile(argv[5]); // file name of the output file

	// an: used to compute thetaW
	float an = 0.0;
	float an2 = 0.0;

	for(i=1; i<nIndiv; ++i){
		an += (1.0 / i);
		an2 += (1.0 / (i*i));
	}


	// define the boundaries of the surveyed windows over a sequence
	std::vector <float> minBin;
	std::vector <float> maxBin;
	window(WIDTH, STEP, minBin, maxBin);
	const size_t nBins(minBin.size());

	// prepare the output file
	write_header(outputFile, nBins);

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
		float pi_avg_bins[nBins];
		std::string haplotypes[nIndiv];*/
		std::vector <unsigned> valide_dataset; // size = nReplicate
		std::vector <unsigned> nSNPtotal; // size = nReplicate
		std::vector < std::vector <float> > pi_avg_bins; // size = [nBins][nReplicate]
		std::vector < std::vector <float> > pi_std_bins; // size = [nBins][nReplicate]
		std::vector < std::vector <float> > thetaW_bins; // size = [nBins][nReplicate]
		std::vector < std::vector <float> > tajimaD_bins; // size = [nBins][nReplicate]

		std::vector <std::string> haplotypes; // size = nIndiv haplotypes with all SNPs
		std::vector <std::string> bin_haplotypes; // size = nIndiv haplotypes with SNPs within bin
		std::vector <float> positions; // size = nSNPs

		// statistics: arrays of length nBins containing statistics averaged over replicates
		float pi_avg_bin[nBins]; // pi_avg 
		float pi_std_bin[nBins]; // pi_std
		float thetaW_bin[nBins]; // thetaW
		float tajimaD_bin[nBins]; // tajD
		
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
					pi_avg_bins.clear();
					pi_std_bins.clear();
					thetaW_bins.clear();
					tajimaD_bins.clear();

					// prepare pi_avg_bins: pi_avg_bins[bin][replicate]
					for(i=0; i<nBins; ++i){
						pi_avg_bins.push_back(tmp_vector_float);
						pi_avg_bin[i] = 0.0;
						pi_std_bins.push_back(tmp_vector_float);
						pi_std_bin[i] = 0.0;
						thetaW_bins.push_back(tmp_vector_float);
						thetaW_bin[i] = 0.0;
						tajimaD_bins.push_back(tmp_vector_float);
						tajimaD_bin[i] = 0.0;
						for(j=0; j<nReplicate; ++j){
							pi_avg_bins[i].push_back(0.0);
							pi_std_bins[i].push_back(0.0);
							thetaW_bins[i].push_back(0.0);
							tajimaD_bins[i].push_back(0.0);
						}
					}

//					std::cout << "set of param: " << cntSetOfParam << std::endl;
				}

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

//				std::cout << "\treplicate_" << replicateID << "\tnSNP " << nSNPtotal[replicateID] << "\tvalidate" << valide_dataset[replicateID] << std::endl;


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

//						std::cout << "Set of param: "<< cntSetOfParam << "\tReplicate: " << replicateID << "\tBin: " << i << std::endl;
						//printBins(bin_haplotypes);

						pi_avg_bins[i][replicateID] = 0.0;
						pi_std_bins[i][replicateID] = 0.0;
						thetaW_bins[i][replicateID] = 0.0;
						tajimaD_bins[i][replicateID] = 0.0;
						computePi(bin_haplotypes, pi_avg_bins, pi_std_bins, thetaW_bins, tajimaD_bins, nIndiv, i, nCombIndiv, replicateID, an, an2);

					} // end of loop over bins


					// if the last replicate among replicates
					if(replicateID == (nReplicate - 1)){ // block for treatment of nReplicate
						// compute statistics for each bin and for each replicate
//						std::cout << "set of param: " << cntSetOfParam << std::endl;
//						std::cout << "compute average pi per bins over replicates" << std::endl;

						// compute statistics averaged over replicates for each bin in arrays of lengths [nBins]
						for(i=0; i<nBins; ++i){
							pi_avg_bin[i] = mean_piStats(pi_avg_bins[i]); // array of length [nBins] containing pi per bin (averaged over reps)
							pi_std_bin[i] = mean_piStats(pi_std_bins[i]); // array of length [nBins] containing pi_std per bin (averaged over reps)
							thetaW_bin[i] = mean_piStats(thetaW_bins[i]); // array of length [nBins] containing thetaW per bin (averaged over reps)
							tajimaD_bin[i] = mean_piStats(tajimaD_bins[i]); // array of length [nBins] containing tajD per bin (averaged over reps)
						}

						write_results(outputFile, pi_avg_bin, pi_std_bin, thetaW_bin, tajimaD_bin, nBins);

						// display average "pi (thetaW), " for each bin
//						for(i=0; i<nBins; ++i){
//							std::cout << pi_avg_bin[i] << " (" << thetaW_bin[i] << "), "; 
//						}
//						std::cout << std::endl << std::endl;

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

	for(j=0; j<nIndiv; ++j){
		bin_haplotypes.push_back("");
	}

	for(size_t i(0); i<positions.size(); ++i){
		if(test_bin == 1 && positions[i] > maxBin[binID]){
//			std::cout << "stop A: " << positions[i] << " greater than " << maxBin[binID] << std::endl;
			break;
		}
		if(positions[i] > minBin[binID] && positions[i] <= maxBin[binID]){
			test_bin = 1;
			for(j=0; j<nIndiv; ++j){
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
	std::cout << std::endl;
}


void computePi(std::vector<std::string> const & bin_haplotypes, std::vector< std::vector< float> > & pi_avg_bins, std::vector< std::vector< float> > & pi_std_bins, std::vector< std::vector< float> > & thetaW_bins, std::vector< std::vector< float> > & tajimaD_bins, size_t const nIndiv, size_t const binID, const unsigned nCombIndiv, unsigned int replicateID, float const an, float const an2){
	// takes as entry bin_haplotypes: a vector of vector containing alignments of nIndiv haplotypes for each the nBins bins
	// results stored in different: vector< vector< float> > 
	// compute pi_avg for the bin 'binID' (rep replicateID), recorded in pi_avg_bins[binID][replicateID]
	// compute pi_std for the bin 'binID' (rep replicateID), recorded in pi_std_bins[binID][replicateID]
	// compute thetaW for the bin 'binID' (rep replicateID), recorded in thetaW_bins[binID][replicateID]

	size_t const sizeBin(bin_haplotypes[0].size());
	float tmpPi(0.0);
	float meanPi(0.0);
	std::string allele="";
	float pi_over_SNPs[sizeBin];
	
	// compute thetaW
	thetaW_bins[binID][replicateID] = sizeBin / an;

//	std::cout << "nSNPs = " << sizeBin << std::endl;

	// if no SNP in the bin --> return(pi = 0.0)
	if(sizeBin == 0){
		pi_avg_bins[binID][replicateID] = 0.0;
		pi_std_bins[binID][replicateID] = 0.0;
		thetaW_bins[binID][replicateID] = 0.0;
		tajimaD_bins[binID][replicateID] = 0.0;
	}else{
		for(size_t i(0); i<sizeBin; ++i){ // loop over SNPs
			tmpPi = 0.0;
			pi_over_SNPs[i] = 0.0;

			for(size_t j(0); j<nIndiv; ++j){ // loop over individuals
				allele = bin_haplotypes[j][i];
				pi_over_SNPs[i] += std::stoi(allele);
			} // end of loop over individuals

//			std::cout << "nAllele 1 = " << pi_over_SNPs[i] << std::endl;
	
			tmpPi = (pi_over_SNPs[i] * (nIndiv - pi_over_SNPs[i])) / nCombIndiv; // tmpPi = #1 * #0 / nCombIndiv
			pi_over_SNPs[i] = tmpPi; // array of pi over SNPs --> use in the futur to compute pi_std and thetaW
//			std::cout << "pi = " << tmpPi << std::endl;
			meanPi += tmpPi; // sum of pi values over SNPs
		} // end of loop over SNPs

		// compute std
		float pi_std(0.0);
		for(size_t i(0); i<sizeBin; ++i){ // loop over SNPs to compute std
			pi_std += (pi_over_SNPs[i] - meanPi/sizeBin) * (pi_over_SNPs[i] - meanPi/sizeBin); // sum of ( pi_i - mean_pi )²
		}
		
		pi_avg_bins[binID][replicateID] = meanPi ;
		pi_std_bins[binID][replicateID]	= sqrt(pi_std);

		// Tajima D: tajimaD(unsigned const nIndiv, float const pi, float const thetaW, float const an, float an2, const unsigned nS){
		tajimaD_bins[binID][replicateID] = tajimaD(nIndiv,  pi_avg_bins[binID][replicateID], thetaW_bins[binID][replicateID], an, an2, sizeBin);

//		std::cout << "pi = " << pi_avg_bins[binID][replicateID] << "\tthetaW = " << thetaW_bins[binID][replicateID] << "\ttajD = " << tajimaD_bins[binID][replicateID] << std::endl;
	}
}


float mean_piStats(std::vector<float> const & liste){
	// takes as entry array of values
	// return mean value of a statistics over replicates by rejecting unvalidated datasets (i.e: without SNP)
	float res(0.0);
	size_t cnt(0);

	for(size_t i(0); i<liste.size(); ++i){
		++cnt;
		res += liste[i];
	}
	return(res/cnt);
}


void write_header(std::string const outputFile, size_t const nBins){
	// function to write the output file containing all statistics
	// Pi_std
	std::ofstream outputFlux(outputFile.c_str(), std::ios::out);

	if(outputFlux){
		for(size_t i=0; i<nBins; ++i){
			outputFlux << " stat_piAvg_" << i; 
		}
		
		for(size_t i=0; i<nBins; ++i){
			outputFlux << " stat_piStd_" << i; 
		}
		
		for(size_t i=0; i<nBins; ++i){
			outputFlux << " stat_thetaW_" << i; 
		}
		
		for(size_t i=0; i<nBins; ++i){
			outputFlux << " stat_tajimaD_" << i; 
		}
		outputFlux << std::endl;
	}else{
		std::cerr <<  "ERROR: cannot oppen the file " << outputFile << std::endl;
		exit(0);
	}
}

void write_results(std::string const outputFile, float const pi_avg_bin[], float const pi_std_bin[], float const thetaW_bin[], float const tajimaD_bin[], size_t const nBins){
	// function to write the output file containing all statistics
	// Pi_std
	std::ofstream outputFlux(outputFile.c_str(), std::ios::app);

	if(outputFlux){
		for(size_t i=0; i<nBins; ++i){
			outputFlux << " " << pi_avg_bin[i]; 
		}
		
		for(size_t i=0; i<nBins; ++i){
			outputFlux << " " << pi_std_bin[i]; 
		}
		
		for(size_t i=0; i<nBins; ++i){
			outputFlux << " " << thetaW_bin[i]; 
		}
		
		for(size_t i=0; i<nBins; ++i){
			outputFlux << " " << tajimaD_bin[i]; 
		}
		outputFlux << std::endl;
	}else{
		std::cerr <<  "ERROR: cannot oppen the file " << outputFile << std::endl;
		exit(0);
	}
}


float tajimaD(unsigned const nIndiv, float const pi, float const thetaW, float const an, float an2, const unsigned nS){
	// nInd = number of individuals in the alignment
	// pi = sum(Kxy over SNPs and pairwise comparisons) / number_of_pairwise_comparisons
	// thetaW = nSNPs / an
	// an = sum(1/(1:(nIndiv-1)))
	// an = sum(1/(1:(nIndiv-1))**2)
	// nS = number of segregating sites

	// b1 and b2
	float b1((nIndiv + 1.0) / (3.0 * (nIndiv - 1.0)));
	float b2(2.0 * (nIndiv*nIndiv + nIndiv + 3.0) / (9.0*nIndiv * (nIndiv - 1.0)));

	// c1 and c2
	float c1(b1 - 1.0 / an);
	float c2(b2 - (nIndiv + 2.0) / (an * nIndiv) + an2/(an*an));
	
	// e1 and e2
	float e1(c1/an);
	float e2(c2/(an*an + an2));
	
	// denominateur
	float denominator(e1*nS + e2*nS*(nS - 1.0));
	denominator = sqrt(denominator);

	// tajima D
	return((pi - thetaW) / denominator);

}

