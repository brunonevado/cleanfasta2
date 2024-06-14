//
//  gff.cpp
//  cleanFasta2
//
//  Created by Bruno on 30/10/2023.
//

#include "gff.hpp"
#include "common.h"


void gff::readFile(bool include, const std::string & scaffold, const std::string & feature, bool verbose){
    if(verbose){
        std::cout << "Reading gff file " << this->infile << " ..." << std::flush;
    }
    this->include = include;
    this->scaffold = scaffold;
    this->feature = feature;
    
    if(feature == "all"){
        this->readAll(verbose);
    }
    else{
        this->readFeature(feature, verbose);
    }
    if(verbose){
        std::cout << " done" << std::endl;
    }
}


void gff::readAll(bool verbose){
    std::string cline;
    unsigned int featuresRead = 0;
    unsigned int sitesMatch = 0;
    std::ifstream readFile(this->infile);
    if(!readFile.is_open()){
        throw std::string("ERROR: cant open gff infile?");
    }
    while (getline (readFile, cline)) {
        std::vector <std::string> fields = msplit(cline, "\t");
        if(fields.at(0) != this->scaffold){
            continue;
        }else{
            featuresRead++;
            try{
                for(unsigned int i = std::stoi(fields.at(3)); i <= std::stoi(fields.at(4)); i++){
                    // note the -1 below, so the positions become 0-indexed
                    this->regionsOfInterest[i-1] = true;
                    sitesMatch++;
                }
            }catch(...){
                std::cout << std::endl << cline << std::endl;
                throw std::string("Error reading gff line");
            }
        }
    }
    if(verbose){
        std::cout << " read " << featuresRead << " matching features (all, " << sitesMatch << " sites) from gff file " << this->infile << " ...";
    }
}

void gff::readFeature(std::string feature, bool verbose){
    std::string cline;
    unsigned int featuresRead = 0;
    unsigned int sitesMatch = 0;
    std::ifstream readFile(this->infile);
    if(!readFile.is_open()){
        throw std::string("ERROR: cant open gff infile?");
    }
    while (getline (readFile, cline)) {
        std::vector <std::string> fields = msplit(cline, "\t");
        if(fields.at(0) != this->scaffold || fields.at(2) != feature){
            continue;
        }else{
            featuresRead++;
            try{
                for(unsigned int i = std::stoi(fields.at(3)); i <= std::stoi(fields.at(4)); i++){
                // note the -1 below, so the positions become 0-indexed
                this->regionsOfInterest[i-1] = true;
                    sitesMatch++;
                }
            }catch(...){
                std::cout << std::endl << cline << std::endl;
                throw std::string("Error reading gff line");
            }
        }
    }
    if(verbose){
        std::cout << " read " << featuresRead << " matching features (" << feature << ", " << sitesMatch << " sites) from gff file " << this->infile << " ...";
    }
}


std::vector <unsigned int> gff::filterSitesByGff (std::vector <unsigned int> & idxs, bool verbose){
    std::vector <unsigned int> res;
    if(this->include){
        for (auto &i : idxs ){
            if(this->regionsOfInterest.count(i) !=0){
                res.push_back(i);
            }
        }
    }
    else{
        for (auto &i : idxs ){
            if(this->regionsOfInterest.count(i) ==0){
                res.push_back(i);
            }
        }
    }
    return res;
}



/*
std::vector <int> gff::parseCDS( int start, int end, std::string strand, int frame){
    std::vector < int> res;

    
}
*/
