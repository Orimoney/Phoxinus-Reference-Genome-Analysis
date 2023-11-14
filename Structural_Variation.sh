##Software requirements include:
                  #-Minimap2
                  #-Syri
                  #-Plotsr

refgenome="Phoxinus_Haplotype_1.fa"
qrygenome="Phoxinus_Haplotype_2.fa"
SAM="Phoxinus_SV.sam"

# Run Minimap2 with the following options:
# -ax: specify the alignment preset to use (in this case, "asm5" recommended for mapping two closely related fasta files)
# -t: specify the number of threads to use
# -o: specify the output file name and path
# --eqx: Output =/X CIGAR operators for sequence match/mismatch, important for compartibility with Syri
# --secondary=no: do not output secondary alignments


minimap2 -ax asm5 --eqx -t 60 $refgenome $qrygenome -o $SAM --secondary=no

#Run Syri wit the following options
# -k : Keeps intermediate files
# -F : Input format with S indicating Sam format i.e tyhe output of the previous minimap alignment

syri -c $SAM -r $refgenome -q $qrygenome -k -F S

#Plot Structural variants with Plotsr
# --sr output of Syri run
# --genomes list of both haplotypes and their paths
# -H height of png output
# -W Width of png output
# -o Output png files
# --chrord preferred order of Chr from a chromosome list

plotsr --sr syri.out  --genomes genomes.txt -H 10 -W 4 -o Phoxinus_Structural_Variants.png --chrord Chrorder.txt
