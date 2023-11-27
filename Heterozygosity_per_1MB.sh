#!/bin/bash
##Software requirements include:
                  #-Minimap2
                  #-Samtools/1.10
                  #-Sambamba
                  #-Bedtools/2.29.2
                  #-ANGSD


################################################ Map Hifi fastq file to the assembly , Sort and remove duplicates ######################################################################
######## Map raw hifi reads to reference genome#################################
#Input files
Reference="Phoxinus_haplotype1.fa" #Reference genome in fasta format
Hifi_fq="Phoxinus_Hifi.fq" #Hifi raw reads in fastq format
Sam_file="Phoxinus.sam" # alignment output file
Bam_file="Phoxinus.bam" # Binary alignment output file
Sorted_Bam_file="Phoxinus.sorted.bam"
Sorted_deduped_Bam_file="Phoxinus.sorted.rm.bam"


# Run Minimap2 with the following options:
# -ax: specify the alignment preset to use (in this case, "asm20" recommended for mapping pacbio long reads to reference genome)
# -t: specify the number of threads to use
# -o: specify the output file name and path
# --secondary=no: do not output secondary alignments
minimap2 -ax asm20 -t 32 ${Reference} ${Hifi_fq} -o ${Sam_file} --secondary=no

# Convert .sam file to .bam file to save storage space using samtools view with the follwoing parameters
# -b : to indicate output is in .bam format
# -q : mapping quality score filter

samtools view -b -q 20  ${Sam_file} > ${Bam_file}

#remove .sam files tro save space
rm ${Sam_file}

#Sort .bam files
# -T: Temporary directory and prefix
# -o: Output file
samtools sort ${Bam_file} -T ../tmp/phoxinus -o ${Sorted_Bam_file}


#Remove duplicates from sorted .bam files
# -r: flag to remove duplicate instead of just flagging them
# --overflow-list-size: increasing this reduces size of temporary files created, very important with a limited cluster
sambamba markdup -r -t 40 ${Sorted_Bam_file} ${Sorted_deduped_Bam_file} --overflow-list-size 800000  --tmpdir=.../tmp


###################################################### Runs of Heterozygosity ###################################
# files
Reference="Phoxinus_haplotype1.fa"
Sorted_deduped_Bam_file="Phoxinus.sorted.rm.bam"
OUTFILE=Phoxinus
SCAFFOLDLIST=${OUTFILE}_1MB_scaffold.lengths.txt # see below for how to generate list
WINDOWSLIST=${OUTFILE}_1MB_windows.txt

# variables
THREADS=10 # number of threads
MEANDEPTH=30 # has to be an integer for mindepth/maxdepth calc below.
	let MINDEPTH=MEANDEPTH/3  # 1/3x average coverage
	let MAXDEPTH=MEANDEPTH*2  # 2x average coverage
MBQ=20  # minimum base quality filter
MAPQ=30  # minimum map quality filter

# Index your fasta file and .bam file
samtools faidx ${Reference}

samtools index ${Sorted_deduped_Bam_file}


# Generate lists of scaffold lenths from the indexed genome fasta files
awk '$2 > 1000000 {print $1"\t"$2}' ${OUTFILE}_haplotype1.fai > ${OUTFILE}_1MB_scaffold.lengths.txt

# Generate 1MB sliding windows from previously generated scaffold lists
bedtools makewindows -g ${OUTFILE}_1MB_scaffold.lengths.txt -w 1000000 | awk '$3 ~ "000000" {print$1":"$2"-"$3}' > ${OUTFILE}__1MB_windows.txt

# Get working scaffold based on array number
echo "Task id is $SGE_TASK_ID"  #To ensure parallel running of multiple chromosomes in a sun grid engine

NUM=$(printf %02d ${SGE_TASK_ID})

# generate windows list for each scaffold
CHR=$(head -n ${NUM} ${SCAFFOLDLIST} | tail -n 1)
CHR=$(echo ${CHR} | awk -F " " '{ print $1 }')
CHR=$(echo ${CHR} | awk -F " " '{$0=$0":"}{print}')
echo ${CHR}

grep ${CHR} ${WINDOWSLIST} > /${OUTFILE}_${NUM}_windows.txt

################################################################
# Run angsd
################################################################
# .saf files are put into temp directory (scratch) to speed up processing and reduce clutter in output directory. See http://www.popgen.dk/angsd/index.php/Input#Genotype_Likelihood_Files for explanation of options.
angsd -GL 2 -setMinDepth ${MINDEPTH} -setMaxDepth ${MAXDEPTH} -minmapq ${MAPQ} -minq ${MBQ} -uniqueonly 1 -remove_bads 1 -only_proper_pairs 1 -docounts 1 -i ${BAMFILE} -ref ${REF} -P ${THREADS} -out ${TEMPDIR}/${OUTFILE}_${NUM} -doSaf 1 -anc ${Reference} -rf ${OUTFILE}_${NUM}_windows.txt -baq 2

# Run realSFS
################################################################
while read -r line; do realSFS -r $line ${TEMPDIR}/${OUTFILE}_${NUM}.saf.idx  -fold 1 -P ${THREADS} -tole 1e-8 2> log >> ${OUTFILE}_1MB_combined_sfs_${NUM}.txt; done < ${Reference}_${NUM}_windows.txt
# -P = number of threads

################################################################
# Calculate heterozygosity (ends up being the value in column 4)
################################################################
awk '{print $1,$2,$3=$1+$2,$4=$2/$3}'${OUTFILE}_1MB_combined_sfs_${NUM}.txt >> ${OUTFILE}_1MB_ROH_${NUM}_summary.het
