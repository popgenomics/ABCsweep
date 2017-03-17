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
	void initialise(); // initialise and clear all attributes
	void ajouteReplicat(unsigned nIndiv); // add one more replicate
	void ajouteSNP(unsigned nSNP); // add the objerved number of SNP to m_nSNP
	void ajoutePosition(std::vector<float> listePos); // add the list of positions between 0 and 1 
	void ajouteHaplotype(std::vector<std::string> listeSeq); // add the list of sequences 
	void setOfParamID(int setID);
	void afficherContenu();

	// attributes:
	private:
	int m_setID; // ID of the replicated set of parameters
	int m_nReplicate; // number of stored replicates of a set of parameters (i.e: with nSNP>0)
	std::vector<unsigned> m_nIndiv; // if 3 replicates : m_nIndiv = [nIndiv_rep1, nIndiv_rep1, nIndiv_rep1]
	std::vector<unsigned> m_nSNP; // if 3 replicates : m_nSNP = [nSNP_rep1, nSNP_rep1, nSNP_rep1]
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
//	bins.printBins();

	// read the msms outputfile
	std::ifstream fifo(msmsFile.c_str());
	if(fifo){ // if the file msmsFile exists
		std::string ligne;
		std::string mot;
		unsigned i(0);
		unsigned nDataset(0); // count the number of simulated datasets over the whole msmsFile
		unsigned nSNPs(0); // number of SNPs for a given dataset
//		int test(-1); // =-1 at the 'segsites' line. setted to +1 at the 'positions' line. haplotypes are only readen if test==1
		unsigned cntIndiv(0); // count the number of haplotypes after the 'positions' line for a given dataset
		std::vector<std::string> haplotypes; // contains nIndiv haplotypes
		Alignment data;
		unsigned replicateID(0);

		while(std::getline(fifo, ligne)){ // read the msmsFile
			// 'segsites' line
			if(ligne[0]=='s'){ // if the line starts by a 's', then we expect: 'segsites: int'
				nDataset++;
				replicateID = nDataset%nReplicate; // the ID of the dataset over nReplicate replicates
				if(replicateID == 1){
				// if first dataset of the replicated combination of parameters
					data.initialise();
				}
//				test = -1;
				std::istringstream iss(ligne);
				i = 0;
				// read over the line containing 'segsites' to get the number of SNPs
				while( std::getline( iss, mot, ' ') ){
					if( i==1 ){ // records the 2nd item of the line as the number of SNPs
						nSNPs = std::stoul(mot);
					}
					i++;
				}
				if(nSNPs > 0){
					data.ajouteReplicat(nIndiv);
					data.ajouteSNP(nSNPs);
				}
				continue; // stop reading the 'segites' line
			} // end of: "if the line contains the string: 'segsites'

		
			// 'positions' line	
			if(ligne[0]=='p'){ // if the line starts by a 'p', then we expect: 'positions: ...'
//				test = 1;
				cntIndiv = 0;
				i = 0; // count the number of elements in the line
				haplotypes.clear(); // clear the vector containing nIndiv sequences
				
				std::istringstream iss(ligne);
				std::vector<float> positions; // vector containing the positions
				// read over the line containing 'positions' to get the positions
				while( std::getline( iss, mot, ' ') ){
					if( i>0 ){ // if mot is not the first item of the line:
						positions.push_back(std::stof(mot));
					}
					i++;
				} // stop reading the 'positions' line
				data.ajoutePosition(positions);
			continue;
			} // end of: "if the line contains the string: 'positions'

			/* If the line is not empty, and doesn't start by '/', nor 'segsites', nor 'positions'*/
			//if(test == 1 && ligne!="" && ligne[0]!='/' && ligne[0]!='s' && ligne[0]!='p'){ 
			if(ligne!="" && ligne[0]!='/' && ligne[0]!='s' && ligne[0]!='p'){ 
				if(cntIndiv < (nIndiv+1)){
					cntIndiv++;
					haplotypes.push_back(ligne); // add the ligne containing a haplotype to the vector <haplotypes>
				}
				if(cntIndiv == nIndiv){
					for(i=0; i<nIndiv; i++){
						data.ajouteHaplotype(haplotypes); // put vector <haplotypes> to Alignment data
					}
/*					std::cout << "end of treatment of dataset " << nDataset <<
					" from replicate " << replicateID << std::endl ;*/
					if(replicateID == 0){
						data.setOfParamID(nDataset/nReplicate);
						data.afficherContenu();
					}
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
Alignment::Alignment() : m_nReplicate(0), m_setID(-1)
{
}

void Alignment::initialise(){
	m_setID = -1;
	m_nReplicate = -1;
	m_nIndiv.clear();
	m_nSNP.clear();
	m_dataset.clear();
	m_position.clear();
}

void Alignment::ajouteReplicat(unsigned nIndiv){ 
	m_nReplicate++ ;
	std::vector <std::string> tmp_seq;
	std::vector <float> tmp_pos;
	m_nIndiv.push_back(nIndiv); // vector of nReplicate values of nIndiv
	m_nSNP.push_back(0); // vector of nReplicate values of nSNPs
	m_dataset.push_back(tmp_seq); // rReplicate vectors, each containing nIndiv sequences
	m_position.push_back(tmp_pos); // nReplicate vectors, each containing nSNP positions
}

void Alignment::ajouteSNP(unsigned nSNP){
	m_nSNP[m_nReplicate] = nSNP;
}


void Alignment::ajoutePosition(std::vector<float> listePos){
	m_position[m_nReplicate] = listePos;	
}


void Alignment::ajouteHaplotype(std::vector<std::string> listeSeq){
	m_dataset[m_nReplicate] = listeSeq;
}


void Alignment::setOfParamID(int setID){
	m_setID = setID;
}


void Alignment::afficherContenu(){
	if( m_nReplicate < 0 ){
		std::cerr << "Nothing to display for combination of parameters " << std::endl;
		exit(0);
	}
	unsigned i(0);
	unsigned j(0);
	std::cout << "Replicated set of parameters #_" << m_setID << " --> " << m_nReplicate + 1 << " replicates" << std::endl;
	for(i=0; i<= m_nReplicate; i++){ // 1) loop over replicates
		std::cout << "replicate " << i << " : " << m_nSNP[i] << " SNPs" << std::endl;
		std::cout << "Positions: ";
		for(j=0; j<m_nSNP[i]; j++){ // 2) loop over SNPs
			std::cout << m_position[i][j] << " ";
		}
		std::cout << std::endl;

		for(j=0; j<m_nIndiv[i]; j++){
			std::cout << m_dataset[i][j] << std::endl;
		}
		std::cout << std::endl;

	}
}



