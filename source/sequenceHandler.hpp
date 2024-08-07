//
//  sequenceHandler.hpp
//  cleanFasta2
//
//  Created by Bruno on 18/07/2023.
//

#ifndef sequenceHandler_hpp
#define sequenceHandler_hpp

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <ctime>
#include "gff.hpp"
#include <sstream>
//using namespace std;


const unsigned int maxseqs = 500;
const unsigned int maxlen = 80000000;

struct Alignment {
    char seqs[maxseqs][maxlen]; // this will hold the sequence data
    int missing[maxseqs][maxlen]; // and this will hold 0s and 1s for missing/present
};


class sequenceHandler {
    
private:
    Alignment* alignedData;
    std::vector < std::string > SeqNames; //sample names as they appear in input file
    unsigned int actualNbp = 0;
    unsigned int actualNSeqs = 0;
    std::vector < float > fractionSequencedPerSite; //for each site, the proportion of samples sequenced
    std::string mainProgramVersion = "N/A";
    std::string inputFile = "N/A";
    bool hasGff = false;

    
    // this function will define what are missing chars (N,n,-)
    bool isMissing(char in);
   
    // returns true if base is one of: a, c, g, t
    bool isHomozygous(char in);

    // will read sequence names from file and place them on vector referenced in call
    void readSampleNames(std::vector < std::string > & SampleNames, std::string infile, bool verbose = false);
    
    // will fill up the vector SampleIndexes with the indexes of samples matching the names in SampleNames
    // leaves checking whether all samples were found to the user (compare size of 2 vectors after running)
    void findIndexesOfSamples(std::vector < std::string> & SampleNames, std::vector <unsigned int> & SampleIndexes , bool verbose = false);
    
    // will fill up the vector referenced with indexes of sites with missing data below threshold
    // must have called calcSamplesPerSite
    // can also use windows here
    void filterSitesByMissing (std::vector <unsigned int> & idxs, float minSequenced, unsigned int start, unsigned int stop, bool verbose = false);
    
    // will write to outfile the sequences in iSeqs and the sites in iSites in fasta format
    void genericFastaWriter(std::string outfile, std::vector <unsigned int> &iSeqs, std::vector <unsigned int> &iSites, bool perFeature = false, bool verbose = false);
    
    // will write to outfile the sequences in iSeqs and the sites in iSites in vcf-like format
    void genericVcfWriter(std::string outfile, std::vector <unsigned int> &iSeqs, std::vector <unsigned int> &iSites, bool verbose = false);
    
    // will calculate and write to outfile per sample missingness & nucleotide diversity
    void genericStatsWriter(std::string outfile, std::vector <unsigned int> &iSeqs, std::vector <unsigned int> &iSites, bool perFeature = false, bool verbose = false);

    //
    void writeVcfLine(std::vector <unsigned int> &iSeqs, unsigned int site, std::ofstream &outfile);
    
    // string with current time
    std::string getCurrentTime(){time_t now = time(0); return std::string(ctime(&now));}

    // fills in res vector positions index and index+1 with bases seen in individual idxInd site idxSite
    // (converts from iupac if needed; assumes diploids)
    void fromIUPAC(unsigned int idxInd, unsigned int idxSite, std::vector <char> & res, unsigned int index);
    
    // fills in res vector with (unique) alleles seen in vector alleles (ignores Ns and gaps)
    void getAllAllelesInSite(std::vector <char> & alleles, std::vector <char> & res);
    
    // get genotype for vcf - first allele in vector uniqueAlleles is the ref, remaining alt(s)
    // allAlleles contains all bases seen in site, and will output genotype for individual in idxSeq and idxSeq+1 (diploids)
    std::string getGenotypeForVcf(std::vector <char> uniqueAlleles, std::vector <char> & allAlleles, unsigned int idxSeq);
    
    // check if the last codon refernced in the array contains a stop codon in any sample
    // the array is a 2-vector containing the indexes on the alignedData, and whether they need to be complemented
    int removeLastCodonIfStop(std::vector < std::vector <int> > & CDSmatrix, std::vector < unsigned int > idxsSequencesToConsider, bool verbose = false);
    
    // returns complementary base
    char getComplementary(unsigned int idxInd, unsigned int idxSite);
    
public:
    sequenceHandler( std::string Pversion = "N/A" ) { this->alignedData = new Alignment; this->mainProgramVersion = Pversion;};
    void readFasta (std::string infile, bool verbose = false);
    unsigned int getNumberSeqs(){return this->actualNSeqs;}
    unsigned int getNumberSites(){return this->actualNbp;}
    
    // will find indexes of samples of interest (either all, or only those defined in infile
    void getIndexesOfSamplesOfInterest( std::vector <unsigned int> & SampleIndexes, std::string infile , bool verbose = false);
    
    // will fill up fractionSequencedPerSite (object's own vector), considering only sequences in indexes defined in iSeqsToConsider
    void calcSamplesPerSite(std::vector <unsigned int> & iSeqsToConsider, bool verbose = false);
    
    // will write either windows or full. For window output will add "windowStats.windowEnd.fas" to outfile name
    // will only consider samples in idxsSequencesToConsider and sites with >= minSequenced
    // will also consider the "per feature" option to output stats per feature in gff
    void genericOutput(std::string outfile, std::vector <unsigned int> iSeqs, unsigned int windowSize, float minSequenced, std::string format = "fasta", bool perGffFeature = false, bool verbose = false);
    
    // will "attach" the a gff object
    void addGff(gff * g){this->myGff = g;};
    
    // will read gff
    void readGff(bool include, const std::string & scaffold, const std::string & feature = "all", bool verbose = false){this->hasGff = true; this->myGff->readFile(include, scaffold, feature, verbose);};

    // get CDS indexes et al
    void returnCDS2mergeFromGff( std::vector < std::vector <int> > & CDSmatrix, std::vector < unsigned int > idxsSequencesToConsider, bool verbose = false);
    
    //write merged CDSs
    void writeCDS(std::string outfile, std::vector < std::vector <int> > & CDSmatrix, std::vector < unsigned int > idxsSequencesToConsider, bool verbose = false);

    
    gff * myGff;

    void test();
};



#endif /* sequenceHandler_hpp */

