##Software requirements include:
                  #-BUSCO
                  #-Diamond


#Input files from BRAKER3 protein annotation

GFF_File="Phoxinus_Braker3.gff3"

Protein_sequences="Phoxinus_braker3_sequences.aa"

Hap1="Phoxinus_Haplotype1.fa"

Trembl_fasta="Trembl.fasta" #Downloaded TrEMBL  database Fasta file


#####################################################BUSCO and AGAT sumamry #########################################################################

#BUSCO run with the following Input
        # -m : mode , to be run in "protein" module
        # -l : lineage , the actinopterygii_odb10 lineage was used, specific to fish genomes
        # -o : Output prefix

busco -m protein -i ${Protein_sequences} -o Phoxinus -l actinopterygii_odb10

#AGAT run to summarize statistics of proteins predicted by BRAKER3

agat_sp_statistics.pl -gff ${GFF_File} -g $Hap1 -p -o Hap1_final_Summary.txt


################################ DIAMOND DATABASE HOMOLOGY SEARCH #######################################################################

#Diamond used for fucntional annotation by blasting our predicted proteins against three different databases i.e TrEMBL, Swissprot and PDBAA, using TrEMBL as an example

diamond makedb --in $Trembl_fasta --db TREMBl_Database/Trembl.dmnd #Downloaded trembl fasta file used as input and a database created

diamond blastp --query ${Protein_sequences} --db TREMBl_Database/Trembl.dmnd --out Phoxinus_trembl_blast.txt --outfmt 6 --evalue 0.00001 --threads 30 --tmpdir ~/tmp

#Count number of unique hits from blast results

cut -f 1 Phoxinus_trembl_blast.txt | sort -u | wc -l > Tremblblast_count.txt

#Extract Gene found in each blast
cut -f 1 Phoxinus_trembl_blast.txt | sort | uniq -c |sort -nr -k 1 > Trembl_genes_list.txt
