#!/bin/sh

# 2nd script

# i) It runs bowtie2 to align paired fastq files against the ananas genome, which produces a sam file

# ii) It runs samtools to produce a bam file from the sam file

# iii) It deletes the sam file (to save space)

# iv) It sorts the bam file

# v) It deletes the unsorted bam file (to save space)

# Overall change - new sorted-deduplicated-minimum-data-quality bam files are created

#names='Brom_Garden1 Brom_Garden4 Brom_Garden9 Brom_Garden13 Brom_Garden20 Brom_Garden21 
#Brom_Garden25 Brom_Garden36 Brom_Garden39 Brom_Garden50 Brom_Garden70 Brom_Garden74 Brom_Garden90'

#names='Atre Rspi_B1596 Vitat Wlin_B758 Pmac_B742 Tion_B84 Tjunc Taus_B1293'

names='Trau_B92 Brom_Garden4_dsamp Brom_Garden70_dsamp Blax_B294 Ldod_B127 Pmac_B742 Rspi_B1596 Wlin_B758'

for name in $names
do

echo $name

echo $(date)

# i) align Tillandsia to Ananas

# --very-sensitive-local
# --sensitive-local - this one
# --fast-local

bowtie2 --sensitive-local -p 12 -x ../pineapple/ananas_index \
-1 ../filtering_2022/${name}_no_org_r1.fastq.gz \
-2 ../filtering_2022/${name}_no_org_r2.fastq.gz -S ${name}.aligned.sam

# ii) convert to bam

samtools view -Shb -@ 12 -o ${name}.aligned.bam ${name}.aligned.sam

rm ${name}.aligned.sam

# iii) sort

samtools sort -@ 12 -o ${name}.aligned.sorted.bam ${name}.aligned.bam

rm ${name}.aligned.bam

# iv) deduplicate

samtools rmdup ${name}.aligned.sorted.bam ${name}.aligned.sorted.nodup.bam

rm ${name}.aligned.sorted.bam

# remove low data quality - not needed

#samtools view -b -q 10 -@ 12 ${name}.aligned.sorted.nodup.bam -o ${name}.aligned.sorted.nodup.MQ10.bam
    
#rm ${name}.aligned.sorted.nodup.bam

# vi) index the BAM files (so that they can be viewed with IGV) - done later

# samtools index -b ${name}.sens_aligned.sorted.nodup.bam

done

###########################################





