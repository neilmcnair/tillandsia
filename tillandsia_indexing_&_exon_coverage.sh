#!/bin/sh

# 3rd script

# to be run in the folder with bam files, and ananas files in sister folder

# Overall change - index files are created, coverage files are created, average coverages are calculated

#names='Atre Brom_Garden4 Brom_Garden9 Brom_Garden13 Brom_Garden20 
#Brom_Garden21 Brom_Garden25 Brom_Garden36 Brom_Garden39 Brom_Garden50 Brom_Garden70 
#Brom_Garden74 Brom_Garden90 Blax Ldod Pmac Rspi Taus Tcur Tdem Tdis Tdiv 
#Tion Tjunc Trau Vitat Wlin'

names='Atre'

for name in $names
do

echo $name

echo $(date)


# i) Index the BAM files (so that they can be viewed with IGV)


samtools index -b ${name}.aligned.sorted.nodup.bam


# ii) Calculate the coverage for each exon
    

    bedtools coverage -sorted -bed -mean \
    -a $REPEX/alignment/Ananas_alignment/pineapple.20150427.exons.sorted_v3.gff3 \
    -b $REPEX/coverage_check/${name}.aligned.sorted.nodup.bam \
    > $REPEX/coverage_check/${name}.aligned.sorted.nodup.exon_coverage.no_genome.bed


# iii) Calculate the average coverage for each read


# proper weighting
awk '{weight=$5-$4+1} {sum+=$10*weight} {leng+=weight} END {printf "%f\n", sum/leng}' \
${name}.aligned.sorted.nodup.exon_coverage.no_genome.bed > ${name}.aligned.sorted.nodup.exon_coverage.no_genome.txt


echo $(date)

done

