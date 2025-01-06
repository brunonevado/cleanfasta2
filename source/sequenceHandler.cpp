//
//  sequenceHandler.cpp
//  cleanFasta2
//
//  Created by Bruno on 18/07/2023.
//

#include "sequenceHandler.hpp"

void sequenceHandler::readFasta (std::string infile, bool verbose){
    // Note this function only reads desinterleaved fasta files.
    int cind=0;
    std::string cline;
    this->inputFile = infile;
    if (verbose){
        std::cout << "Reading file " << infile << "..." << std::flush ;
    }
    
    std::ifstream readFile(infile);
    if(!readFile.is_open()){
        throw std::string("ERROR: cant open infile?");
    }
    while (getline (readFile, cline)) {
        if(cind >= maxseqs){
            throw std::string("ERROR: too many sequences, need to recompile program");
        }
        if( cline.at(0) == '>' ){
            this->SeqNames.push_back(cline.substr(1));
            cind++;
        }
        else{
            if(this->actualNbp == 0){
                if(cline.size() > maxlen){
                    throw std::string("ERROR: sequences too long, need to recompile program");
                }
                this->actualNbp = (unsigned int)cline.size();
            }
            for(unsigned int i = 0; i < cline.size(); i++){
                this->alignedData->seqs[cind-1][i] = cline.at(i);
                this->alignedData->missing[cind-1][i] = (cline.at(i) == 'N' || cline.at(i) == '-' || cline.at(i) == 'n' ) ? 0 : 1;
            }
        }
    }
    this->actualNSeqs = cind;
    if(cind >= maxseqs){
        throw std::string("ERROR: too many sequences, need to recompile program");
    }
    readFile.close();
    if (verbose){
        std::cout << " done, read " << this->actualNSeqs << " individuals, " << this->actualNbp << " bps" << std::endl;
    }
}

void sequenceHandler::calcSamplesPerSite(std::vector <unsigned int> & iSeqsToConsider, bool verbose){
    if(verbose){
        std::cout << "Calculating percentage samples sequenced per site considering " << iSeqsToConsider.size() << " samples..." << std::flush;
    }
    this->fractionSequencedPerSite.resize(this->actualNbp);
    // Calculate number of samples present for each position
    for(unsigned int i = 0; i < this->actualNbp; i++){
        for(auto &j : iSeqsToConsider){
            this->fractionSequencedPerSite.at(i) += this->alignedData->missing[j][i];
        }
        this->fractionSequencedPerSite.at(i) /= (float) iSeqsToConsider.size();
    }
    if(verbose){
        std::cout << " done" << std::endl;
    }
}


void sequenceHandler::readSampleNames(std::vector < std::string > & SampleNames, std::string infile, bool verbose){
    
    std::ifstream readFile(infile);
    if(!readFile.is_open()){
        throw std::string("ERROR: cant open infile " + infile + " with names of samples of interest");
    }
    std::string cline;
    while (getline (readFile, cline)) {
        if(cline.size() > 0) {
            SampleNames.push_back(cline);
        }
    }
    if(verbose) {
        std::cout << "Read " << SampleNames.size() << " sequence names to retrieve" << std::endl;
    }
    
    
}




/*
 void sequenceHandler::writeFasta(string outfile, bool filterMissing, float thr, bool verbose){
 
 ofstream writeFile(outfile);
 if(!writeFile.is_open()){
 throw string ("ERROR: cant open file for output?");
 }
 
 // output with missing data
 if(filterMissing){
 if(thr < 0 || thr > 1){
 throw string ("ERROR: maxMissing incorrectely defined?");
 }
 if(verbose){
 cout << "Writing to file " << outfile << " sites with at least " << thr*100 << "% samples sequenced..." << flush;
 }
 for(unsigned int i = 0; i < this->actualNSeqs; i++){
 writeFile << this->SeqNames.at(i) << endl;
 for(unsigned int j = 0; j < this->actualNbp; j++){
 if(!filterMissing || this->fractionSequencedPerSite.at(j) >= thr){
 writeFile << this->alignedData->seqs[i][j];
 }
 }
 writeFile << endl;
 }
 if(verbose){
 cout << " done" << endl;
 }
 }
 writeFile.close();
 }
 */

void sequenceHandler::filterSitesByMissing (std::vector <unsigned int> & idxs, float minSequenced, unsigned int start, unsigned int stop, bool verbose){
    if(verbose){
        std::cout << "Filtering sites for missingness..." << std::flush;
    }
    for(unsigned int i = start; i < stop; i++){
        if(this->fractionSequencedPerSite.at(i) >= minSequenced){
            idxs.push_back(i);
        }
    }
    if(verbose){
        std::cout << " done, found " << idxs.size() << " sites to output" << std::endl;
    }
}

void sequenceHandler::findIndexesOfSamples(std::vector < std::string> & SampleNames, std::vector <unsigned int> & SampleIndexes , bool verbose ){
    if(verbose){
        std::cout << "Finding position of " << SampleNames.size() << " samples in data... "<< std::flush;
    }
    for( auto &name : SampleNames ){
        bool alreadyFound = false;
        for( unsigned int i = 0; i < this->getNumberSeqs(); i++){
            if( name == this->SeqNames.at(i) ){
                if(alreadyFound){
                    std::cout << std::endl;
                    throw std::string("ERROR: sample found multiple times – " + name);
                }else{
                    SampleIndexes.push_back(i);
                    alreadyFound = true;
                }
            }
        }
    }
    if(verbose){
        std::cout << " done" << std::endl;
    }
}


void sequenceHandler::getIndexesOfSamplesOfInterest( std::vector <unsigned int> & SampleIndexes, std::string infile , bool verbose ){
    
    if(infile == ""){ // use all samples if -samples not defined
        SampleIndexes.resize(this->getNumberSeqs());
        iota(SampleIndexes.begin(), SampleIndexes.end(), 0);
    }
    else{ // get sample names, then get their indexes in the matrix. Finally, check all samples were found (could change this last behaviour to be more permissive if needed)
        std::vector < std::string > sampleNames;
        this->readSampleNames(sampleNames, infile, verbose);
        this->findIndexesOfSamples(sampleNames, SampleIndexes, verbose);
        
        if(sampleNames.size() != SampleIndexes.size()){
            if(this->strictNames == false){
                std::cerr << "WARNING: not all sequences found in file.";
            }
            else{
                throw std::string("Error: not all sequences found in file. Check input and options or rerun with -strictNames 0.");
            }
        }
    }
}


void sequenceHandler::genericFastaWriter(std::string outfile, std::vector <unsigned int> &iSeqs, std::vector <unsigned int> &iSites, bool perFeature, bool verbose){
    // write single fasta file
    if(!perFeature){
        std::ofstream writeFile(outfile);
        if(!writeFile.is_open()){
            throw std::string ("ERROR: cant open file for output?");
        }
        if(verbose){
            std::cout << "Writing fasta format to file " << outfile << " : " << iSeqs.size() << " sequences, " << iSites.size() << " sites..." << std::flush;
        }
        for( auto &i : iSeqs){
            writeFile << ">" << this->SeqNames.at(i) << std::endl;
            for(auto &j : iSites){
                writeFile << this->alignedData->seqs[i][j];
            }
            writeFile << std::endl;
        }
            writeFile.close();
    }
    // write fasta file per feature
    else{
        for(int iFeature = 0; iFeature < myGff->getNumberOfFeatures(); iFeature++){
            std::vector <unsigned int> SitesWithinFeature;
            int startFeature = myGff->getFeatureStart(iFeature);
            int stopFeature = myGff->getFeatureStop(iFeature);
            for(auto &idxSite: iSites){
                if(idxSite < startFeature){
                    continue;
                }
                else if(idxSite > stopFeature){
                    break;
                }
                else{
                    SitesWithinFeature.push_back(idxSite);
                }

            }
            // note here we call this same function but with the perFeature set to false, so it outputs the current feature
            if(SitesWithinFeature.size() > 0){
                std::string outfileFeature = outfile;
                outfileFeature.append(".feature");
                outfileFeature.append(std::to_string(iFeature));
                outfileFeature.append(".fas");
                this->genericFastaWriter(outfileFeature, iSeqs, SitesWithinFeature,false, verbose);
                SitesWithinFeature.clear();
            }
            
        }
    }
    if(verbose){
        std::cout << " done" << std::endl;
    }
}


void sequenceHandler::genericOutput(std::string outfile, std::vector <unsigned int> iSeqs, unsigned int windowSize, float minSequenced, std::string format, bool perGffFeature, bool verbose){
    
    // output - all
    if(windowSize <= 0){
        std::vector < unsigned int > sitesToOutput;
        this->filterSitesByMissing(sitesToOutput, minSequenced, 0, this->getNumberSites(), verbose);
        
        if(this->hasGff){
            std::vector <unsigned int> filteredSites = this->myGff->filterSitesByGff(sitesToOutput);
            sitesToOutput.clear();
            sitesToOutput = filteredSites;
        }
        
        if(format == "fasta") {
            this->genericFastaWriter(outfile, iSeqs, sitesToOutput, perGffFeature, verbose);}
        else if (format == "vcf"){
            this->genericVcfWriter(outfile, iSeqs, sitesToOutput, verbose);
        }
        else if (format == "stats"){
            this->genericStatsWriter(outfile, iSeqs, sitesToOutput, perGffFeature, verbose);
        }
    }
    // output - windowed
    else{
        unsigned int start=0;
        while(start < this->getNumberSites()){
            std::vector < unsigned int > sitesToOutput;
            unsigned int stop = start + windowSize;
            if(stop > this->getNumberSites()){
                stop = this->getNumberSites();
            }
            std::string outfileWindow = outfile + "." + std::to_string(start+1) + "_" + std::to_string(stop);
            this->filterSitesByMissing(sitesToOutput, minSequenced, start, stop, verbose);
            
            if(this->hasGff){
                std::vector <unsigned int> filteredSites = this->myGff->filterSitesByGff(sitesToOutput);
                sitesToOutput.clear();
                sitesToOutput = filteredSites;
            }
            if(format == "fasta") {
                outfileWindow.append(".fas");
                this->genericFastaWriter(outfileWindow, iSeqs, sitesToOutput,perGffFeature, verbose);
            }
            else if (format == "vcf"){
                outfileWindow.append(".vcf");
                this->genericVcfWriter(outfileWindow, iSeqs, sitesToOutput, verbose);
            }
            else if (format == "stats"){
                outfileWindow.append(".stats");
                this->genericStatsWriter(outfileWindow, iSeqs, sitesToOutput, perGffFeature, verbose);
            }
            else{
                throw std::string("Unknown format specified (this is a coding issue?): " + format);
            }
            start+= windowSize;
        }
    }
}


void sequenceHandler::genericVcfWriter(std::string outfile, std::vector <unsigned int> &iSeqs, std::vector <unsigned int> &iSites, bool verbose){
    std::ofstream writeFile(outfile);
    if(!writeFile.is_open()){
        throw std::string ("ERROR: cant open file for output?");
    }
    if(verbose){
        std::cout << "Writing vcf format to file " << outfile << " : " << iSeqs.size() << " sequences, " << iSites.size() << " sites..." << std::flush;
    }
    
    // VCF header
    writeFile << "##fileformat=VCFv4.x" << std::endl << "##File created with cleanfasta2 v" << this->mainProgramVersion << " from input fasta file " << this->inputFile << " on " << this->getCurrentTime();
    writeFile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    writeFile << "##INFO=<ID=PS,Number=1,Type=Float,Description=\"Proportion sequenced\">" << std::endl;
    // eventually would like to add the call to the program, similar to: ##bcftools_viewCommand=view Lup2342.unflt.bcf; Date=Fri Jul 28 12:47:25 2023
    writeFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    
    for( auto &i : iSeqs){
        writeFile << "\t" << this->SeqNames.at(i);
    }
    writeFile << std::endl;
    
    // Now process each site at a time
    for(auto &i : iSites){
        try{
            this->writeVcfLine(iSeqs, i, writeFile);
        }catch(...){
            throw std::string("Unable to process site " + std::to_string(i+1) + " into vcf format");
        }
    }
    writeFile << std::endl;
    if(verbose){
        std::cout << " done" << std::endl;
    }
    writeFile.close();
}

void sequenceHandler::writeVcfLine(std::vector <unsigned int> &iSeqs, unsigned int site, std::ofstream &outfile){
    
    // create vector with all alleles seen in position (from IUPAC if needed)
    std::vector <char> allAlleles (iSeqs.size()*2, '0');
    unsigned int posInAlleleVector = 0;
    for(auto &i: iSeqs){
        this->fromIUPAC(i, site, allAlleles, posInAlleleVector );
        posInAlleleVector += 2;
    }
    // get unique alleles seen in position
    std::vector <char> uniqueAlleles;
    this->getAllAllelesInSite(allAlleles, uniqueAlleles);
    
    // write out beginning of line
    outfile << this->inputFile << "\t" << site+1 << "\t.\t" << uniqueAlleles.at(0) << "\t";
    // write out alternative alleles seen
    if(uniqueAlleles.size() == 1){
        outfile << ".";
    }else{
        outfile << uniqueAlleles.at(1);
        for(unsigned int i = 2; i < uniqueAlleles.size(); i++){
            outfile << "," << uniqueAlleles.at(i);
        }
    }
    // write out info section + GT section
    outfile <<  "\tPASS\tPS=" << this->fractionSequencedPerSite.at(site) << "\tGT";
    // write out genotypes
    for(unsigned int i = 0; i < allAlleles.size(); i +=2){
        outfile << "\t" << this->getGenotypeForVcf(uniqueAlleles, allAlleles, i);
    }
    
    outfile << std::endl;
}

void sequenceHandler::fromIUPAC(unsigned int idxInd, unsigned int idxSite, std::vector <char> & res, unsigned int index){
    switch(this->alignedData->seqs[idxInd][idxSite]){
        case 'A':
            res[index]='A';
            res[index+1]='A';
            break;
        case 'C':
            res[index]='C';
            res[index+1]='C';
            break;
        case 'G':
            res[index]='G';
            res[index+1]='G';
            break;
        case 'T':
            res[index]='T';
            res[index+1]='T';
            break;
        case 'Y':
            res[index]='C';
            res[index+1]='T';
            break;
        case 'K':
            res[index]='G';
            res[index+1]='T';
            break;
        case 'S':
            res[index]='C';
            res[index+1]='G';
            break;
        case 'W':
            res[index]='A';
            res[index+1]='T';
            break;
        case 'R':
            res[index]='A';
            res[index+1]='G';
            break;
        case 'M':
            res[index]='A';
            res[index+1]='C';
            break;
        case 'N':
            res[index]='N';
            res[index+1]='N';
            break;
        default:
            std::cerr << "Unexpected IUPAC code: " << this->alignedData->seqs[idxInd][idxSite] << std::endl;
            throw std::string("Received unexpected IUPAC code.");
    }
}

void sequenceHandler::getAllAllelesInSite(std::vector <char> & alleles, std::vector <char> & res){
    
    for(auto &i: alleles){
        if(i == 'N' || i == 'n'){
            continue;
        }
        bool seen = false;
        for(auto &j: res){
            if(j == i){
                seen = true;
                continue;
            }
        }
        if(!seen){
            res.push_back(i);
        }
    }
}

std::string sequenceHandler::getGenotypeForVcf(std::vector <char> uniqueAlleles, std::vector <char> & allAlleles, unsigned int idxSeq){
    // missing - note I assume if missing both bases are missing
    if(allAlleles.at(idxSeq) == 'N' || allAlleles.at(idxSeq) == 'n'){
        return "./.";
    }
    else{
        std::vector <unsigned int> genotype;
        for(unsigned int i = 0; i < uniqueAlleles.size(); i++){
            if(allAlleles.at(idxSeq) == uniqueAlleles.at(i) || allAlleles.at(idxSeq+1) == uniqueAlleles.at(i)){
                genotype.push_back(i);
            }
        }
        std::string res = std::to_string(genotype.at(0));
        res += "/";
        res += std::to_string(genotype.at(genotype.size()-1));
        return res;
        
    }
    return "";
}

void sequenceHandler::genericStatsWriter(std::string outfile, std::vector <unsigned int> &iSeqs, std::vector <unsigned int> &iSites, bool perGffFeature, bool verbose){
    if(verbose){
        std::cout << "Calculating and outputting stats for " << iSeqs.size() << " samples, " << iSites.size() << " sites...";
    }
    std::ofstream writeFile(outfile);
    if(!writeFile.is_open()){
        throw std::string ("ERROR: cant open file for output?");
    }
    // stats for everything in gff
    if(!perGffFeature){
        for(auto &idxInd: iSeqs){
            unsigned int seqSites = 0;
            unsigned int iupacSites = 0;
            for(auto &idxSite: iSites){
                if(!isMissing (this->alignedData->seqs[idxInd][idxSite] )){
                    seqSites++;
                    if(!isHomozygous(this->alignedData->seqs[idxInd][idxSite] )){
                        iupacSites++;
                    }
                }
            }
        writeFile << this->SeqNames.at(idxInd) << "\t" << 1-  float(seqSites) / float(iSites.size()) << "\t" << float(iupacSites) / float(seqSites) <<  std::endl;
        }
    }
    // write stats per feature
    else{
        for(int iFeature = 0; iFeature < myGff->getNumberOfFeatures(); iFeature++){
            for(auto &idxInd: iSeqs){
                unsigned int seqSites = 0;
                unsigned int iupacSites = 0;
                for(auto &idxSite: iSites){
                    if (idxSite < myGff->getFeatureStart(iFeature)){
                        continue;
                    }
                    else if (idxSite > myGff->getFeatureStop(iFeature)){
                        break;
                    }
                    if(!isMissing (this->alignedData->seqs[idxInd][idxSite] )){
                        seqSites++;
                        if(!isHomozygous(this->alignedData->seqs[idxInd][idxSite] )){
                            iupacSites++;
                        }
                    }
                }
            writeFile << this->SeqNames.at(idxInd) << "\t" << 1-  float(seqSites) / float(myGff->getFeatureStop(iFeature) - myGff->getFeatureStart(iFeature) +1 )  << "\t" << float(iupacSites) / float(seqSites) <<  std::endl;
            }
        }
    }
    writeFile.close();

    if(verbose){
        std::cout << " done" << std::endl;
    }
}


bool sequenceHandler::isMissing(char in){
    if( in == 'N' || in == 'n' || in == '-'){
        return true;
    }else{
        return false;
    }
}

bool sequenceHandler::isHomozygous(char in){
    if( in == 'A' || in == 'T'  || in == 'G' || in == 'C'  || in == 'a' || in == 't'  || in == 'g' || in == 'c' ){
        return true;
    }else{
        return false;
    }
}


void sequenceHandler::test(){
    //this->myGff->readFile(true,"scaffold_1","Gypsy_LTR_retrotransposon",true);
    this->myGff->test();
}

// writes just the cds. still need to complement if needed
void sequenceHandler::writeCDS(std::string outfile, std::vector < std::vector <int> > & CDSmatrix, std::vector < unsigned int > idxsSequencesToConsider, bool verbose ){
    std::ofstream writeFile(outfile);
    if(!writeFile.is_open()){
        throw std::string ("ERROR: cant open file for output?");
    }
    if(verbose){
        std::cout << "Writing fasta format to file " << outfile  << std::flush;
    }
    for(auto & i: idxsSequencesToConsider){
        writeFile << ">" << this->SeqNames.at(i) << std::endl;
        for(int idxCDS = 0; idxCDS < CDSmatrix.at(0).size(); idxCDS++){
            if(CDSmatrix.at(1).at(idxCDS) == 1){
                writeFile << getComplementary(i, CDSmatrix.at(0).at(idxCDS));
            }
            else{
                writeFile << this->alignedData->seqs[i][CDSmatrix.at(0).at(idxCDS)];
            }
        }
        writeFile << std::endl;
    }
    if(verbose){
        std::cout << " done" << std::endl;
    }
    writeFile.close();
    
}

// gets CDS positions
void sequenceHandler::returnCDS2mergeFromGff( std::vector < std::vector <int> > & CDSmatrix, std::vector < unsigned int > idxsSequencesToConsider, bool verbose){
    // This will return the index (0-based) of sites in CDS, in the correct order for each CDS (i.e., from end to start if CDS is in reverse strand)
    // note we skip the last codon if incomplete, and also the first couple bases if out of frame
    // will also fill in a second vector with 0 for not need to complement, and 1 if needed to complement
    // finally, will exclude last codon if it is a stop codon.
    std::vector < int> res;
    std::string feature = "CDS";
    std::string cline;
    unsigned int featuresRead = 0;
    unsigned int sitesMatch = 0;
    unsigned int stopCodonsRemoved = 0;
    std::ifstream readFile(this->myGff->getInfile());
    if(!readFile.is_open()){
        throw std::string("ERROR: cant open gff infile?");
    }
    while (getline (readFile, cline)) {
        std::vector <std::string> fields = msplit(cline, "\t");
        if(fields.at(0) != this->myGff->getScaffold() || fields.at(2) != feature){
            continue;
        }else{
            featuresRead++;
            try{
                //fields.at(3) = start; std::stoi(fields.at(4)) = end; fields.at(6) = strand; fields.at(7) = frame
                // gff values are 1-based
                if( fields.at(6) == "+" ){
                    // CDS in forward strand
                    int start = std::stoi(fields.at(3)) + std::stoi(fields.at(7)) - 1;
                    int end = std::stoi(fields.at(4)) - 1;
                    for(int i = start; i < end; i += 3){
                        //last codon incomplete - just skip
                        if( i + 3 > end){
                            break;
                            
                        }
                        else{
                            CDSmatrix.at(0).push_back(i);
                            CDSmatrix.at(0).push_back(i+1);
                            CDSmatrix.at(0).push_back(i+2);
                            CDSmatrix.at(1).push_back(0);
                            CDSmatrix.at(1).push_back(0);
                            CDSmatrix.at(1).push_back(0);
                            sitesMatch += 3;
                        }
                    }
                }
                else{
                    // CDS in reverse strand
                    int start = std::stoi(fields.at(4)) - std::stoi(fields.at(7)) - 1;
                    int end = std::stoi(fields.at(3)) - 1;
                    for(int i = start; i > end; i -= 3){
                   //last codon incomplete - just skip
                        if( i - 3 < end){
                            break;
                        }
                        else{
                            CDSmatrix.at(0).push_back(i);
                            CDSmatrix.at(0).push_back(i-1);
                            CDSmatrix.at(0).push_back(i-2);
                            CDSmatrix.at(1).push_back(1);
                            CDSmatrix.at(1).push_back(1);
                            CDSmatrix.at(1).push_back(1);
                            sitesMatch += 3;
                        }
                    }
                }
            }
        catch(...){
                std::cout << std::endl << cline << std::endl;
                throw std::string("Error reading gff line");
            }
            // here we have read a CDS. Before continuing we check for stop codons and remove if present
            stopCodonsRemoved += removeLastCodonIfStop( CDSmatrix, idxsSequencesToConsider, verbose);

        }
    }
    if(verbose){
        std::cout << "Merging CDS: read " << featuresRead << " matching features (" << feature << ", " << sitesMatch << " sites) from gff file " << this->myGff->getInfile() << ", removed " << stopCodonsRemoved << " stop codons..." << std::endl;
    }
}

// gets base complementary to base in position idxSite of individual idxInd
char sequenceHandler::getComplementary(unsigned int idxInd, unsigned int idxSite){
    switch(this->alignedData->seqs[idxInd][idxSite]){
        case 'A':
            return('T');
            break;
        case 'C':
            return('G');
            break;
        case 'G':
            return('C');
            break;
        case 'T':
            return('A');
            break;
        case 'Y':
            return('R');
            break;
        case 'K':
            return('M');
            break;
        case 'S':
            return('S');
            break;
        case 'W':
            return('W');
            break;
        case 'R':
            return('Y');
            break;
        case 'M':
            return('K');
            break;
        case 'N':
            return('N');
            break;
        default:
            std::cerr << "Unexpected IUPAC code: " << this->alignedData->seqs[idxInd][idxSite] << std::endl;
            throw std::string("Received unexpected IUPAC code.");
    }
    
}

// checks if any of the samples of interest has a stop codon on the last codon added to CDSmatrix
// if so removes it, and returns 1; else returns 0;
int sequenceHandler::removeLastCodonIfStop(std::vector < std::vector <int> > & CDSmatrix, std::vector < unsigned int > idxsSequencesToConsider, bool verbose){
    // need to loop over each individual to consider
    for (auto &idxInd : idxsSequencesToConsider){
        std::string lastCodon;
        // get last codon
        for(int i = CDSmatrix.at(0).size() - 3; i < CDSmatrix.at(0).size(); i++){
            if( CDSmatrix.at(1).at(i) == 0){
                lastCodon.push_back(this->alignedData->seqs[idxInd][CDSmatrix.at(0).at(i)]);
            }
            else{
                lastCodon.push_back( getComplementary(idxInd, CDSmatrix.at(0).at(i)));
            }
        }
        // now we remove the last 3 bases and break if it matches stop codon. Note:using standard genetic code here
        if (lastCodon == "TAG" || lastCodon == "TAA" || lastCodon == "TGA" ){
            std::cout << " found stop codon: " << lastCodon << std::endl;
            // Dont know what would be fastest here
            CDSmatrix.at(0).pop_back();
            CDSmatrix.at(0).pop_back();
            CDSmatrix.at(0).pop_back();
            CDSmatrix.at(1).pop_back();
            CDSmatrix.at(1).pop_back();
            CDSmatrix.at(1).pop_back();
            return 1;
            
        }
    }
    
    return 0;
}


