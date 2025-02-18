//
//  main.cpp
// cleanFasta2
//  Created by Bruno on 17/07/2023.
//

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <numeric>

#include "args2.h"
#include "gff.hpp"
#include "sequenceHandler.hpp"


std::string Pversion = "v2.060125";

void help(){
    std::cout << "###################\n  cleanFasta2 " << Pversion << " \n###################" << std::endl;
    std::cout << "Will read fasta file and write new fasta file filtering sites and sequences according to different criteria" << std::endl;
    std::cout << "Usage example: cleanFasta2 -infile infile.fas -outfile outfile.fas -maxMissing 50 [-verbose 1] [ -winsize 1000000 ] [ -samples samples.txt ] [-format fasta/vcf/stats] [-gff annot.gff -scaffold scaf1 -include 0 -feature repeat_region -perFeature 0] [-mergeCDS 0] [-strictNames 1]" << std::endl;
    std::cout << "-infile: input fasta file, should be desinterleaved!" << std::endl;
    std::cout << "-outfile: output fasta file. If using windows or gff, will be used as prefix." << std::endl;
    std::cout << "-maxMissing: only output sites with missingness % at most maxMissing (considering only some samples if -samples set)." << std::endl;
    std::cout << "-winsize: output in windows of defined size (unset or set to < 0 to output single window)." << std::endl;
    std::cout << "-samples: text file containing samples to consider (1 sample per line - if not used outputs all samples)." << std::endl;
    std::cout << "-verbose: some info to screen if set to != 0 (default = 0, i.e. no screen output)." << std::endl ;
    std::cout << "-format: output format, either fasta (the default), vcf or stats." << std::endl << std::endl;
    std::cout << "-strictNames: if set to 1 (the default), will exit program with an error if not all samples in -samples file are found in infile. If set to 0, will just raise a warning and continue." << std::endl << std::endl;
    
    std::cout << "GFF options:" << std::endl;
    std::cout << "-gff: GFF annotation file with regions to include or exclude." << std::endl;
    std::cout << "-scaffold: name of scaffold in GFF file that matches input file." << std::endl;
    std::cout << "-include: whether to include (set to 1) or exclude (set to 0) the features in GFF." << std::endl;
    std::cout << "-feature: name of features to include or exclude (by default will include/exclude all features in gff)." << std::endl;
    std::cout << "-mergeCDS: use in conjunction with gff and CDS regions. Will complement/reverse and place in frame each CDS, and output all CDS (in each window if selected). Ignores missingness filter." << std::endl;
    std::cout << "-perFeature: use in conjunction with gff. With stats output will write missingness and nucleotide diversity per sample and per feature (same order as in gff). With fasta output will write 1 fasta file per feature (named outfile.featureX.fas, where X is position of feature in GFF file)." << std::endl << std::endl;

    std::cout << "NOTES: 1. Only reads desinterleaved fasta format" << std::endl;
    std::cout << "       2. If -samples used, all samples must be present in fasta file" << std::endl;
    std::cout << "       3. Program will read entire matrix to memory, and at start requests enough memory for 2 matrixes (1 int, 1 char) of 500 samples x 80 MB sequences." << std::endl;
    std::cout << "          These values can be changed in sequenceHandler.h (requires recompilation). " << std::endl;
    std::cout << "       4. Ignores gaps (actually, probably crashes with gaps under some settings)." << std::endl;
    std::cout << "       5. Output format 'stats' will write out missingness and nucleotide diversity per sample (after excluding sites with missing data above treshold, and by windows if -winsize used ." << std::endl;

    
    
}

/*
 Still to add:
 1. read interleaved file
 . check if input file is aligned
 */

// 220524: added check for when matrix is too large
// 030424: Added a different analysis - merge CDSs. This will ignore all other options for filtering except "samples"
// it exports the alignment from the CDSs present in GFF file, reverse complement them if needed, and put in phase. Then remove
// stop codons if at end, and finally concatenate everything into 1 long CDS-like
// it needs to have the gff and the corresponding scaffold set in call to program
// 010724: added perFeature option, which will output missingness and Pi for each selected sample (after missingness filter) in each selected GFF feature
// 290724: changed checks for too long/too many sequences
// 070824: added option to output per feature fasta file
// 060125: added option not to crash if not all samples are found

int main(int argc, const char * argv[]) {
    
    // read command line options
    margs programOptions;
    try{
        programOptions.getargs(argc, argv, std::vector <std::string> {
            "infile,s,f",
            "outfile,s,f",
            "maxMissing,i,f",
            "samples,s,t",
            "winsize,i,t",
            "verbose,b,t",
            "format,s,t",
            "gff,s,t",
            "scaffold,s,t",
            "include,b,t",
            "feature,s,t",
            "perFeature,b,t",
            "strictNames,b,t",
            "mergeCDS,b,t"});
    }catch(std::string e){ help();std::cerr << std::endl << "Failed reading args: " << e << std::endl;exit(1);}

    std::string infile = programOptions.getString("infile");
    std::string outfile = programOptions.getString("outfile");
    float maxMissing = (float)programOptions.getInt("maxMissing")/100.0;
    float minSequenced = 1 - maxMissing;
    unsigned int winsize = programOptions.isArgDefined("winsize") ? programOptions.getInt("winsize") : 0;
    std::string sampleNamesFile = programOptions.isArgDefined("samples")  ? programOptions.getString("samples") : "";
    std::string outputFormat = programOptions.isArgDefined("format")  ? programOptions.getString("format") : "fasta";
    if(outputFormat != "vcf" && outputFormat != "fasta" && outputFormat != "stats"){
        std::cerr << "Unknown output format requested: " << outputFormat << " (allowed are fasta, vcf and stats)" << std::endl;
        exit(1);
    }
    bool noisy = programOptions.isArgDefined("verbose") ? programOptions.getBool("verbose") : false;
    // gff options & checks
    bool mergeCDS = programOptions.isArgDefined("mergeCDS") ? programOptions.getBool("mergeCDS") : false;
    bool strictNames = programOptions.isArgDefined("strictNames") ? programOptions.getBool("strictNames") : false;
    std::string gffFile = programOptions.isArgDefined("gff")  ? programOptions.getString("gff") : "";
    std::string gffScaffold = programOptions.isArgDefined("scaffold")  ? programOptions.getString("scaffold") : "";
    std::string gffFeature = programOptions.isArgDefined("feature")  ? programOptions.getString("feature") : "all";
    bool gffInclude = programOptions.isArgDefined("include")  ? programOptions.getBool("include") : false;
    if(programOptions.isArgDefined("gff") && (!programOptions.isArgDefined("scaffold") || !programOptions.isArgDefined("include"))){
        std::cerr << "For GFF filtering -scaffold and -include must to be defined" << std::endl;
        exit(1);
    }
    
    if(mergeCDS && gffFile == "" ){
        std::cerr << "For mergeCDS option need to provide a GFF file" << std::endl;
        exit(1);
    }
    

    // option to output stats per feature in gff file (instead of stats for concatenation of feature)
    // if this option is true, output must be GFF - needs check here before proceeding; and window options must not be set; and include must be 1; and output must be stats
    bool gffPerFeature = programOptions.isArgDefined("perFeature")  ? programOptions.getBool("perFeature") : false;
    if(gffPerFeature){
        if(gffFile == ""){
            std::cerr << "For perFeature option need to provide a GFF file" << std::endl;
            exit(1);
        }
        if(winsize != 0){
            std::cerr << "The perFeature option is incompatible with winsize != 0" << std::endl;
            exit(1);
            }
        if(!gffInclude){
            std::cerr << "The perFeature option is incompatible with include = 0" << std::endl;
            exit(1);
        }
        if(outputFormat != "stats" && outputFormat != "fasta"){
            std::cerr << "The perFeature option only available with stats or fasta output formats" << std::endl;
            exit(1);
        }
    }



    // create objects
    sequenceHandler mySeqHandler (Pversion);
    gff agff(gffFile);
    mySeqHandler.addGff(&agff);
    mySeqHandler.setStrictNames(strictNames);
     
    try{
        // read fasta file - desinterleaved only for now!
        mySeqHandler.readFasta(infile, noisy);
        
        // check small size, indicative of interleaved format
        if (mySeqHandler.getNumberSites() < 500){
            std::cout << "<cleanfasta2> WARNING: sequences are short (< 500 bp), is the input fasta file interleaved?" << std::endl;
        }

        // read gff if defined
        if(gffFile != ""){
            mySeqHandler.readGff(gffInclude,gffScaffold, gffFeature,noisy);
        }
        
        
        // get indexes of samples of interest into idxsSequencesToConsider (either all samples or just those in sample infile
        std::vector < unsigned int > idxsSequencesToConsider;
        mySeqHandler.getIndexesOfSamplesOfInterest(idxsSequencesToConsider, sampleNamesFile, noisy);
        
        // calculate per site missingness for sequences of interest
        mySeqHandler.calcSamplesPerSite(idxsSequencesToConsider, noisy);
        

        // if mergeCDS option, here just get the indexes of sites to export
        if( mergeCDS == true){
            std::vector< std::vector <int> > CDSs;
            CDSs.resize(2);
            mySeqHandler.returnCDS2mergeFromGff(CDSs,idxsSequencesToConsider, noisy);
            mySeqHandler.writeCDS(outfile, CDSs,idxsSequencesToConsider, noisy);
        }
        else{
            // output
            mySeqHandler.genericOutput(outfile, idxsSequencesToConsider, winsize, minSequenced, outputFormat, gffPerFeature, noisy);
        }
    }
    catch(std::string e){
        std::cerr << e << std::endl;
        exit(1);
    }
    
    
    
    
    
    return 0;
}

