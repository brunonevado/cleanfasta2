# cleanfasta2
Re-implementation of cleanfasta to clean fasta files according to different criteria.

This version is much faster but requires a lot more memory.  
It does not use Bio++ libraries so should be easier to compile.  


Will read fasta file and write new fasta file filtering sites and sequences according to different criteria  
Usage example: cleanFasta2 -infile infile.fas -outfile outfile.fas -maxMissing 50 [-verbose 1] [ -winsize 1000000 ] [ -samples samples.txt ] [-format fasta/vcf/stats] [-gff annot.gff -scaffold scaf1 -include 0 -feature repeat_region -perFeature 0] [-mergeCDS 0]  
  
-infile: input fasta file, should be desinterleaved!  
-outfile: output fasta file. If using windows or gff, will be used as prefix.  
-maxMissing: only output sites with missingness % at most maxMissing (considering only some samples if -samples set).  
-winsize: output in windows of defined size (unset or set to < 0 to output single window).  
-samples: text file containing samples to consider (1 sample per line - if not used outputs all samples).  
-verbose: some info to screen if set to != 0 (default = 0, i.e. no screen output).  
-format: output format, either fasta (the default), vcf or stats.  
  
    
GFF options:  
-gff: GFF annotation file with regions to include or exclude.  
-scaffold: name of scaffold in GFF file that matches input file.  
-include: whether to include (set to 1) or exclude (set to 0) the features in GFF.  
-feature: name of features to include or exclude (by default will include/exclude all features in gff).  
-mergeCDS: use in conjunction with gff and CDS regions. Will complement/reverse and place in frame each CDS, and output all CDS (in each window if selected). Ignores missingness filter.  
-perFeature: use in conjunction with gff, will output missingness and nucleotide diversity per sample and per feature (same order as in gff; only available with stats output).  
  
    
NOTES: 1. Only reads desinterleaved fasta format  
       2. If -samples used, all samples must be present in fasta file.  
       3. Program will read entire matrix to memory, and at start requests enough memory for 2 matrixes (1 int, 1 char) of 500 samples x 80 MB sequences.  
          These values can be changed in sequenceHandler.h (requires recompilation).  
       4. Ignores gaps (actually, probably crashes with gaps under some settings).  
       5. Output format 'stats' will write out missingness and nucleotide diversity per sample (after excluding sites with missing data above treshold, and by windows if -winsize used .  

  
  
Installation (wihtin /source – and need c++17):    
  
g++ -O3 -std=c++17 *cpp -o cleanfasta2

