# Phoxinus Reference Genome



## Phoxinus Reference genome analysis scripts
This directory contains scripts used in the analysis of the Phoxinus phoxinus assembly and annotation manuscript, with descriptions of parameters and tools used for each analysis.

### 1. Busco_Diamond_blast.sh
### Description:
This script performs three post-protein-coding gene prediction analyses, including BUSCO completion assessment, summary of annotation statistics and finally fucntional annotation of predicted protein-coding genes using the following tools:
1. BUSCO
2. AGAT
3. DIAMOND

Required inputs for these analyses are the genome assembly fasta file, the selected database for the diamond blast, in this case TrEMBL, which can be downloaded from https://www.uniprot.org, the .gff3 and protein sequences from the BRAKER3 output.

## 2. Gene_Evolution_Analysis.sh
## Description
This script performs gene evolution analysis with ORTHOFINDER and CAFE5, by identifying the longest isoforms of genes in a .gff3 file and then running it through ORTHOFINDER and CAFE5 to identify orthogroups between selected species and infer gene expansion and contraction relative to other species.

## 3. Heterozygosity_per_1MB and Heterozygosity_Estimates.Rscript
## Description
This script calculates heterozygosity estimates per 1MB window over the entire genome. It divides the genome into specified windows and computes heterozygosity metrics for each window and an R script designed to visualise comprehensive heterozygosity estimates across the entire genome.

## 4. Kmer_Genome_Profile
## Description
This script profiles the genome using k-mer analysis, generates k-mer frequency profiles and gives a preliminary insight into the raw data before genome assembly.

## 5. PSMC.sh
## Description
This PSMC (Pairwise Sequentially Markovian Coalescent) script reconstructs historical population sizes from diploid genome sequences. It uses the PSMC model to infer past demographic events.

## 6. Structural_Variation.sh
## Description
This script identifies and analyzes structural variations within the genome. It leverages various tools to detect deletions, duplications, inversions, and translocations. It begins with  aligning both haplotypes with MINIMAP2 , identifying structural variants with SYRI and plotting the synteny and structural variants with PLOTSR. 







