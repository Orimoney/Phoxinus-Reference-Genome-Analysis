##Software requirements include:
                  #-AGAT
                  #-Orthofinder
                  #-CAFE
                  #-Cafeplotter


GFF_File="Phoxinus_Braker3.gff3"
Longest_Isoform_gff="Phoxinus_longest_Isoform.gff3"
Protein_sequences="Phoxinus_braker3_sequences.aa"
Hap1="Phoxinus_Haplotype1.fa"


#AGAT used to extract  longest isoforms used orthology analysis downstream

agat_sp_keep_longest_isoform.pl -gff ${GFF_File}  -o ${Longest_Isoform_gff}

agat_sp_extract_sequences.pl -g ${Longest_Isoform_gff} -f ${Hap1} -p -o Longest_Isoform_sequences.fa


################################################ ORTHOFINDER run for Orthology #############################################
#This should be run inh the directory conmtaining all proteomes involved in the analysis are
# -f : fasta file
# -S : blast programme to use, in this case DIAMOND
# -M : Model , in thsi case multiple sequence alignment

orthofinder  -f ${Longest_Isoform_sequences} -S diamond -t 40 -M msa

#Convert Orthofinder's species tree result to ultrametric(i.e) species tree with branchlength
# -r : to indicate root age, can be gotten from timetree.org
python /home/toriowo/.conda/envs/AGAT/bin/make_ultrametric.py  -r 142  SpeciesTree_rooted.txt

#Edit the Orthogroups.GeneCount.tsv file to tmp.tsv only extracting Orthogroups and numbers of genes present
awk -F'\t' '{print "(null)\t"$0}' Orthogroups.GeneCount.tsv > tmp.tsv

#Change the first "(null)" header to Desc
sed '0,/null/s//Desc/' tmp.tsv > tmp_file.tsv && mv tmp_file.tsv tmp.tsv

#remove the total column from above, without needed to figure out column numbers.
awk -F'\t' '{$NF=""; print $0}' tmp.tsv | rev | sed 's/^\s*//g' | rev | tr ' ' '\t' > mod.tsv

#filter the Orthogroups.GeneCount.tsv file to remove OG that have more than 100 proteins in a particular species
python2.7 cafetutorial_clade_and_size_filter.py -i mod.tsv -s -o cafe.input.tsv

python2.7 ${Clade_filter} -i mod.tsv -s -o cafe.input.tsv

####################################CAFE5 run to infer gene expansion and contraction in the Phoxinus elative  to the LCA #######################

#Run Cafe5 with ultrametric tree and cafe.input.tsv as input
cafe5 --infile cafe.input.tsv --tree SpeciesTree_rooted.txt.ultrametric.tre

#plot gene expansion/contraction tree with cafeplotter with output of CAFE5
cafeplotter -i results/ -o cafeplot/

#After extracting gene lists of contracting and expanding gene families, merge the sequences together from Orthofinder Orthogroup sequences results
cat $(cat Expansion_Gene_List_Phoxi.txt) | awk '/^>g/{print $0; getline; while(/^\w/){print $0; getline}}' > Expansion_sequences.fa

cat $(cat Gene_Contraction_Phoxinus.txt) | awk '/^>g/{print $0; getline; while(/^\w/){print $0; getline}}' > Contraction_sequences.fa

cat $(cat Phoxinus_specificlist.txt) | awk '/^>g/{print $0; getline; while(/^\w/){print $0; getline}}' > Phoxinus_correct_specific.fa
