######Software requirements
              #Samtools v1.8
              #Bcftools v1.10.2
              #PSMC v0.6.5 conda installation possible
              
#load all 3 programmes either via module load or conda activate as the case may be

#This is important to ensure parallelization of jobs downstream
echo "Task id is $SGE_TASK_ID"

#Paths for input files and output directory, use full paths
H1_bam="Phoxinus_Hap1.sorted.rm.bam"
H1_genome="Phoxinus_Haplotype_1.fa" #Reference genome
output_directory="~/PSMC"

#Declare job to run on arrays, each array per chromosome {This is specific for sungrid engine}, so
#If you specify 50 threads when submitting, it will use 50 threads per chromosome, I suggest using around 15 or 20 threads so 5-8 chromosomes can run at once

CHR="Chr${SGE_TASK_ID}"
#Generate consensus fastq from the .bam files using a combination of samtools mpileup, bcftools call and vcfutils.pl vcf2fq
#This run is in parallel, resulting in multiple consensus fastq files according to the number of chromosomes present
# Q: minimum base quality
# q: minimum map quality
# c: bcftyools consensus caller
# d: minimum depth
# D: Maximum depth

samtools mpileup -Q 30 -q 30 -u -v \
-f ${H1_genome} -r ${CHR} ${H1_bam} | \
bcftools call -c | \
vcfutils.pl vcf2fq -d 10 -D 60 -Q 30 > ${output}/HAP1.${CHR}.fq

#Merge output per chromosome fq files for consensus fq
cat ${output}/*.fq > ${output}/$Hap1.consensus.fq

#Transform the consensus sequence into a fasta like format with PSMC's fq2psmcfa
#If you are not interested in bootstrapping, then you can directly use this .psmcfa file as input to the "psmc" downstream
#-q20: minimum base quality
fq2psmcfa -q20 ${output}/${Hap1}.consensus.fq > $output/Hap1.psmcfa

#If you are however interested in bootstrapping, then psmcfa file has to be split with PSMC splitfa, to make sequence subsampling easier
splitfa Hap1.psmcfa > ${output}/split.psmcfa

#Then finally run the main psmc analysis on 100 randomly sampled sequences from the of the split psmcfa
# -b: ensures random sampling

for i in $(seq 1 100)
do
  psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o $output/Hap1_${i}.psmc split.psmcfa

done

#Merge all 100 bootstrap results in readiness for plotting

cat *.psmc > ${output}/Phoxinus_Bootstrap.psmc

#Plot with psmc_plot.pl
#-p .psmc : output file name
#-g : generation interval in years (depends on your species)
#-u : mutation rate (Can be estimated but usually, mutation rates from close relatives works)
psmc_plot.pl -R -u 3.51e-09 -g 3 -p ${output}/Hap1_Plot_bootstrap Phoxinus_Bootstrap.psmc
