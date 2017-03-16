#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

class Bins{
	// object of class Bins contains two vectors 'min' and 'max'
	// corresponding to the coordinates of all explored windows over a sequence
	// with a given width and step
	public:
	size_t nElements;
	std::vector<float> min;
	std::vector<float> max;
};


class Alignment{
	public:
	int nIndiv;
	int nReplicate;
	std::vector<std::string> dataset; // dataset[over replicates][over individuals] -> one sequence per individual
	
};


Bins window(float width=0.1, float step=0.05);

// arg1 = name of the msms outputfile
int main(int argc, char* argv[]){
	const std::string msmsFile(argv[1]); // mmsFile contains the name of msms output file to read
	const int nIndiv = atoi(argv[2]); // number of individuals in the alignment
	const int nCombParam = atoi(argv[3]); // number of combination of parameters
	const int nReplicate = atoi(argv[4]); // number of replicates per combination of parameters: n sim tot = nCombParam * nReplicate

	// define the boundaries of the surveyed windows over a sequence
	Bins bins;
	bins = window();

	// read the msms outputfile
	std::ifstream fifo(msmsFile.c_str());
	if(fifo){ // if the file msmsFile exists
		std::string ligne;
		std::string mot;
		unsigned i(0);
		unsigned nDataset(0); // count the number of simulated datasets over the whole msmsFile
		unsigned nSNPs(0); // number of SNPs for a given dataset
		int test(-1); // =-1 at the 'segsites' line. setted to +1 at the 'positions' line. haplotypes are only readen if test==1
		unsigned cntIndiv(0); // count the number of haplotypes after the 'positions' line for a given dataset
		std::vector<std::string> alignement;

		while(std::getline(fifo, ligne)){ // read the msmsFile

			// 'segsites' line
			if(ligne[0]=='s'){ // if the line starts by a 's', then we expect: 'segsites: int'
				test = -1;
				std::istringstream iss(ligne);
				i = 0;
				// read over the line containing 'segsites' to get the number of SNPs
				while( std::getline( iss, mot, ' ') ){
					if( i==1 ){ // records the 2nd item of the line as the number of SNPs
						nSNPs = std::stoul(mot);
					}
					i++;
				}
				// stop reading the 'segites' line
//				std::cout << ligne << std::endl;
			continue;
			} // end of: "if the line contains the string: 'segsites'

		
			// 'positions' line	
			if(ligne[0]=='p'){ // if the line starts by a 'p', then we expect: 'positions: ...'
				test = 1;
				cntIndiv = 0;
				i = 0; // count the number of elements in the line
				alignement.clear(); // clear the vector containing nIndiv sequences
				
				std::istringstream iss(ligne);
				std::vector<float> positions; // vector containing the positions
				// read over the line containing 'positions' to get the positions
				while( std::getline( iss, mot, ' ') ){
					if( i>0 ){ // if mot is not the first item of the line:
						positions.push_back(std::stof(mot));
					}
					i++;
				} // stop reading the 'positions' line
			continue;
			} // end of: "if the line contains the string: 'positions'


			if(test == 1 && ligne!="" && ligne[0]!='/' && ligne[0]!='s' && ligne[0]!='p'){
				if(cntIndiv < (nIndiv+1)){
					cntIndiv++;
					alignement.push_back(ligne);
				//	std::cout << ligne << std::endl;
				}
				if(cntIndiv == nIndiv){
					nDataset++;
					for(i=0; i<nIndiv; i++){
						std::cout << i << ": " << alignement[i] << std::endl;
					}
					std::cout << "end of treatment of dataset " << nDataset <<
					" from replicate " << nDataset%nReplicate << std::endl << std::endl;
				}
			}
			
		} // end of loop over the msmsFile
	}else{
		std::cerr << "ERROR: cannot oppen the file " << msmsFile << std::endl;
		exit(0);
	}	
	
	// test
	/* // test for the bins
	Bins test;
	test = window();
	size_t i = 0;
	for(i=0; i<test.nElements; i++){
		std::cout << "min = " << test.min[i] << "\tmax = " << test.max[i] << std::endl;
	}*/
	return(0);
}


Bins window(float width, float step){
	/* Function that return an object of class Bins
	with the positions of the 'min' and 'max' boundaries
	of windows over a sequence */
	Bins res;
	size_t nElements = 1;
	float minB = 0;
	float maxB = width;

	res.min.push_back(minB);
	res.max.push_back(maxB);

	while((minB + width) <= 1 && maxB < 1){
		nElements++;
		minB += step;
		maxB += step;
		if(maxB > 1){
			maxB = 1;
		}
	res.min.push_back(minB);
	res.max.push_back(maxB);
	}
	
	res.nElements = nElements;

	return(res);
}

