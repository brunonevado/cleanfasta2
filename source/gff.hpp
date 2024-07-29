//
//  gff.hpp
//  cleanFasta2
//
//  Created by Bruno on 30/10/2023.
//

#ifndef gff_hpp
#define gff_hpp

#include <iostream>
#include <string>
#include <unordered_map>
#include <fstream>
#include <vector>


#include "common.h"


class gff {
private:
    std::string infile;
    std::string scaffold;
    std::unordered_map <unsigned int, bool> regionsOfInterest;
    bool include;  // true if positions in regionsOfInterest are to be exported false if they are to be excluded
    std::string feature;
    void readAll(bool verbose);
    
    std::vector <int> cdsPositions;
    std::vector <bool> cdsPositionsComplement;

    std::vector <int> featureStartPos;
    std::vector <int> featureStopPos;
    int featuresRead = 0;
    
    void readFeature(std::string feature, bool verbose);

public:
    gff( std::string infile) {this->infile = infile;};
    void readFile(bool include, const std::string & scaffold, const std::string & feature = "all", bool verbose = false);
    // this will iterate of sites in idxs, and output only sites that match gff filter
    std::vector <unsigned int> filterSitesByGff (std::vector <unsigned int> & idxs, bool verbose = false);
    void test(){std::cout << "Gff object has " << this->regionsOfInterest.size() << std::endl;}
    // this will fill up vector with sites to output, and whether to complement them
    //void returnCDS2merge( std::vector < std::vector <int> > & CDSmatrix, bool verbose);
    std::string getInfile(){return infile;}
    std::string getScaffold(){return scaffold;}
    int getNumberOfFeatures(){return featuresRead;}
    int getFeatureStart(int i){return featureStartPos.at(i);}
    int getFeatureStop(int i){return featureStopPos.at(i);}
};






#endif /* gff_hpp */
