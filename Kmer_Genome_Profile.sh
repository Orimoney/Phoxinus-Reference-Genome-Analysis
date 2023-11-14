#Software requirements include:
                    #-Jellyfish/2.3.0 for kmer counting
                    #-genomescope2 2.0 for genome profiling

#Input files
Hifi_fq="Phoxinus_Hifi.fq" #Hifi raw reads in fastq format

#output files
Jf_file="Phoxinus_reads.jf" #output of jellyfish count, contains kmer counts

Histo_file="Phoxinus_reads.histo" #Histogram file of kmer count



#Count k-mers with the following command with -C canonical counting method, -m 19 k-mer length, -t 20 threads and -s 10Gb of memory
jellyfish count -C -m 19 -s 10000000 -t 20 ${Hifi_fq} -o ${Jf_file}


#To export the k-mer count into a histogram, using the output of "jellyfish count above" with -t set to 20 threads

jellyfish histo -t 20 ${Jf_file} > ${Histo_file}


#To run genomescope2 on histo file with -k set to 19 which is the k-mer length used for counting above , will create a folder which contains full reports
genomescope.R -i ${Histo_file} -o output -k 19
