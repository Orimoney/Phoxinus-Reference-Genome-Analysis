# Phoxinus Reference Genome



## Phoxinus Reference genome analysis scripts
This directory contains scripts used in the analysis of the Phoxinus phoxinus assembly and annotation manuscript, with descriptions of parameters and tools used for each analysis.

### 1. Busco_Diamond_blast.sh
### Description:
This script performs three post-protein-coding gene prediction analyses, including BUSCO completion assessment, summary of annotation statistics and finally fucntional annotation of predicted protein-coding genes using the following programmes:
- BUSCO (https://busco.ezlab.org/)
- AGAT  (https://github.com/NBISweden/AGAT)
- DIAMOND (https://github.com/bbuchfink/diamond)

Required inputs for these analyses are the genome assembly fasta file, the selected database for the diamond blast, in this case TrEMBL, which can be downloaded from https://www.uniprot.org, the .gff3 file and protein sequences from the BRAKER3 output.

## 2. Gene_Evolution_Analysis.sh
## Description
This script performs gene evolution analysis with ORTHOFINDER and CAFE5, by identifying the longest isoforms of genes in a .gff3 file and then running it through ORTHOFINDER and CAFE5 to identify orthogroups between selected species and infer gene expansion and contraction relative to other species. The required programmes include:
- BUSCO (https://busco.ezlab.org/)
- ORTHOFINDER (https://github.com/davidemms/OrthoFinder)
- CAFE (https://github.com/hahnlab/CAFE)
- Cafeplotter (https://github.com/moshi4/CafePlotter) 

## 3. Heterozygosity_per_1MB and Heterozygosity_Estimates.Rscript
## Description
This script calculates heterozygosity estimates per 1MB window over the entire genome. It divides the genome into specified windows and computes heterozygosity metrics for each window and an R script designed to visualise comprehensive heterozygosity estimates across the entire genome. It begins by calculating genotype likelihoods and site frequency spectrum (SFS) before estimating heterozygosity from the SFS output. The required programmes include:
- Minimap2 (https://github.com/lh3/minimap2)
- Samtools (https://github.com/samtools)
- Sambamba (https://github.com/biod/sambamba)
- Bedtools (https://bedtools.readthedocs.io/en/latest/index.html#)
- ANGSD (https://www.popgen.dk/angsd/index.php/ANGSD)


  _*Note: The script is written for Sun Grid Engine (SGE). Adjustments may be needed for SLURM schedulers._
 
## 4. Kmer_Genome_Profile
## Description
This script profiles the genome using k-mer analysis, generates k-mer frequency profiles and gives a preliminary insight into the raw data before genome assembly. The required programmes include:
- Jellyfish (https://github.com/gmarcais/Jellyfish)
- GenomeScope2 (https://github.com/tbenavi1/genomescope2.0)

## 5. PSMC.sh
## Description
This PSMC (Pairwise Sequentially Markovian Coalescent) script attempts to reconstruct historical population sizes from diploid genome sequences. It uses the PSMC model to infer past demographic events. It takes a reference genome and .bam file as input. It requirees the following programmes:
- Samtools (https://github.com/samtools)
- Bcftools (https://samtools.github.io/bcftools/)
- PSMC (https://github.com/lh3/psmc)


## 6. Structural_Variation.sh
## Description
This script identifies and analyzes structural variations within the genome. It leverages various tools to detect deletions, duplications, inversions, and translocations. It begins with  aligning both haplotypes with MINIMAP2 , identifying structural variants with SYRI and plotting the synteny and structural variants with PLOTSR.
- Minimap2 (https://github.com/lh3/minimap2)
- SYRI (https://github.com/schneebergerlab/syri)
- PlotSR (https://github.com/schneebergerlab/plotsr)

## 7. genome_assembly.sh
## Description
This script provides a comprehensive pipeline for genome assembly using HiFi reads, followed by duplicate removal, scaffolding, manual curation, and polishing. The prerequuisite tools for this pipeline includ:
- HiFiasm (https://github.com/chhylp123/hifiasm)
- Purge Dups (https://github.com/dfguan/purge_dups)
- SALSA2 (https://github.com/marbl/SALSA)
- BWA (https://github.com/lh3/bwa)
- Samtools (https://github.com/samtools)
- Minimap2 (https://github.com/lh3/minimap2)
- Bedtools ((https://bedtools.readthedocs.io/en/latest/index.html#))
- Pairtools (https://github.com/open2c/pairtools)
- Cooler (https://github.com/open2c/cooler)
- PBMM2 (https://github.com/PacificBiosciences/pbmm2)
- Singularity (https://github.com/sylabs/singularity)
- DeepVariant (https://github.com/google/deepvariant)
- Bcftools (https://samtools.github.io/bcftools/)

   *_Parameters such as the number of threads can be adjusted based on available computational resources._

## 8. Mapping_subreads_hifi_hic.sh
## Description
This script is designed to perform comprehensive mapping and quality assessment of subreads, HiFi, and Hi-C sequencing reads to a reference genome assembly. The script automates the alignment process, merging BAM files, and generating mapping statistics.

## 9. Mercury_genome_completion_check.sh
## Description
This script assesses completion and quality of a genome assembly using Mercury. The follwing programmes are requirements:
- FastK (https://github.com/thegenemyers/FASTK)
- MerquryFK (https://github.com/thegenemyers/MERQURY.FK)




