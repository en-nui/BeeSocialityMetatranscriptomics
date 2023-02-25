# BeeSocialityMetatranscriptomics
Author: Chris Robinson  
Bioinformatic pipeline used for data analysis of RNA viruses of 6 species of bee metatranscriptomics. Also including python script for Fst analysis (it's embarrassing, hello 2nd year of grad school).   

# Bioinformatic Pipeline
Objectives: Viral genome reconstruction, create tables for ecological and evolutionary data, extract viral RdRPs  

Inputs: Bee RNA data from Illumina NovaSeq  
Quick walkthru:  
ASSEMBLY AND VIRAL RECOVERY  
QC reads with MutltiQC -> Trim with Trimmomatic -> Align reads to respective genomes -> Discard aligned reads, keep unaligned reads -> use unaligned reads for assembly with rnaviralSPADES. Use CheckV and Cenote-Taker2 to QC assembly and recover viral contigs.

RdRP SCREENING  
Palmscan and Cenote-Taker2 to recover RdRP sequences from contigs (use this as marker gene for RNA virus in sample, analogus to 16S). 

Phylogeny 
Align RdRP containing contigs with MAFFT and trim with TrimAL. Build tree with IQTREE2. 

#Fst Analysis  
Objective: Investigate variation in RNA virus populations associated with bee hosts

Inputs: Variant calling file (VCF) from LoFreq without Indel calling. from aligning all Narnaviridae contigs in each bee host sample to longest Narnaviridae contig. 

Required Software: 
BWA-mem
LoFreq without Indel Calling

Script
fst_calc.py #for calculation of fst  
pairwise_fst.py #for pairwise calculations between all bee samples  
population_fst.py #for pairwise calculations between all sample sites  

