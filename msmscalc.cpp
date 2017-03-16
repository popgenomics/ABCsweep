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
	//methods:
	public:
	void window(float width, float step);
	void printBins();

	// attributes
	private:
	size_t m_nElements;
	std::vector<float> m_min;
	std::vector<float> m_max;
};


class Alignment{
	// methods:
	public:
	Alignment();
	void ajouteReplicat(); 

	// attributes:
	private:
	unsigned m_nReplicate;
	std::vector<unsigned> m_nIndiv;
	std::vector<unsigned> m_nSNP;
	std::vector< std::vector<std::string> > m_dataset; // dataset[over replicates][over individuals] -> one sequence per individual
	std::vector< std::vector<float> > m_position; // dataset[over replicates][over SNPs] -> one position per SNP
};


//Bins window(float width=0.1, float step=0.05);

// arg1 = name of the msms outputfile
int main(int argc, char* argv[]){
	const std::string msmsFile(argv[1]); // mmsFile contains the name of msms output file to read
	const unsigned nIndiv = atoi(argv[2]); // number of individuals in the alignment
	const unsigned nCombParam = atoi(argv[3]); // number of combination of parameters
	const unsigned nReplicate = atoi(argv[4]); // number of replicates per combination of parameters: n sim tot = nCombParam * nReplicate

	// define the boundaries of the surveyed windows over a sequence
	Bins bins;
	bins.window(0.1, 0.05);
	bins.printBins();

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
		std::vector<std::string> alignement; // contains nIndiv haplotypes

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


			if(test == 1 && ligne!="" && ligne[0]!='/' && ligne[0]!='s' && ligne[0]!='p'){ /* If the line is not
				empty, and doesn't start by '/', nor 'segsites', nor 'positions'*/
				if(cntIndiv < (nIndiv+1)){
					cntIndiv++;
					alignement.push_back(ligne);
				//	std::cout << ligne << std::endl;
				}
				if(cntIndiv == nIndiv){
					nDataset++;
					for(i=0; i<nIndiv; i++){
					//	std::cout << i << ": " << alignement[i] << std::endl;
					}
					std::cout << "end of treatment of dataset " << nDataset <<
					" from replicate " << nDataset%nReplicate << std::endl ;
				}
			}
			
		} // end of loop over the msmsFile
	}else{
		std::cerr << "ERROR: cannot oppen the file " << msmsFile << std::endl;
		exit(0);
	}	
	
	return(0);
}


// methods for class 'Bins'
void Bins::window(float width, float step){
	/* Function that return an object of class Bins
	with the positions of the 'min' and 'max' boundaries
	of windows over a sequence */
	m_nElements = 1;
	float minB = 0.0;
	float maxB = width;

	m_min.push_back(minB);
	m_max.push_back(maxB);

	while((minB + width) <= 1 && maxB < 1){
		m_nElements++;
		minB += step;
		maxB += step;
		if(maxB > 1){
			maxB = 1;
		}
	m_min.push_back(minB);
	m_max.push_back(maxB);
	}
}

void Bins::printBins(){
	size_t i(0);
	for(i=0; i<m_nElements; i++){
		std::cout << "bin #" << i << "\tmin = " << m_min[i] << "\tmax = " << m_max[i] << std::endl;
	}
}


// methods for class 'Alignment'
// constructeur
Alignment::Alignment() : m_nReplicate(0)
{
}

void Alignment::ajouteReplicat(){ 
	std::vector <std::string> tmp_seq;
	std::vector <float> tmp_pos;
	m_nIndiv.push_back(0); // vector of nReplicate values of nIndiv
	m_nSNP.push_back(0); // vector of nReplicate values of nSNPs
	m_dataset.push_back(tmp_seq); // rReplicate vectors, each containing nIndiv sequences
	m_position.push_back(tmp_pos); // nReplicate vectors, each containing nSNP positions
}

