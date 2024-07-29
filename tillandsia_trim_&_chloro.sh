#!/bin/sh

# 1st script

# This script will be run to remove the chloroplast and mitochondrial reads from our raw data (which is fastq.gz)

# i) It runs trimmomatic to make all the reads the same length.

# ii) It runs bowtie2 to align paired fastq files against the ananas chloroplast/mitochondrion genome, which produces a sam file.

# iii) It runs samtools and bedtools to get the output back into fastq.gz format.

# iv) It deletes the intermediate files (to save space).

# Overall change - should produce a pair of fastq.gz files that are smaller than the input, because organelles 
# have been removed and all reads are fixed to be the same length

names='Brom_Garden20 Brom_Garden25 Brom_Garden70'

#names='Blax_B294 Ldod_B127 Pmac_B742 Rspi_B1596 Wlin_B758 Brom_Garden1 Taus_B1293 Tcur_B1597 Tdem_B1373 Tdiv_B1594
# Tjunc Trau_B92 Brom_Garden4 Brom_Garden9 Brom_Garden13 Brom_Garden21 Tion_B84 Brom_Garden90'

for name in $names
do

echo $name

echo $(date)

# trim

java -jar /usr/local/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 12 \
$READDATA/${name}_trim1.fastq.gz $READDATA/${name}_trim2.fastq.gz ${name}_forward_paired.fq.gz \
${name}_forward_unpaired.fq.gz ${name}_reverse_paired.fq.gz ${name}_reverse_unpaired.fq.gz CROP:125 HEADCROP:10 MINLEN:115

# align (using bowtie2)

bowtie2 --very-sensitive-local -p 12 -x ../alignment/Ananas_chloro_mito/indices/Ananas_comosus_chloro_mito_index \
-1 ${name}_forward_paired.fq.gz -2 ${name}_reverse_paired.fq.gz -S ${name}.sam

# filter out mapped reads (using samtools) (so we are left with reads that didn't align to either chloroplast or to mitochondrion)

samtools view -Shb -@ 12 -f 12 -o ${name}.no_organelle.bam ${name}.sam

rm ${name}.sam

# deduplicate

samtools rmdup ${name}.no_organelle.bam ${name}.no_org.nodup.bam

rm ${name}.no_organelle.bam

# filter out mapq < 10 - unnecessary, also lots of reads have mapq = 0, indicating multiple alignment

#samtools view -b -q 10 -@ 12 ${name}.no_org.nodup.bam -o ${name}.no_org.nodup.MQ10.bam
    
#rm ${name}.no_org.nodup.bam

# sort (so that bedtools will work)

samtools sort -n -@ 12 -o ${name}.no_org.nodup.sorted.bam ${name}.no_org.nodup.bam

rm ${name}.no_org.nodup.bam

# convert from bam to fastq (using bedtools) - produces lots of warnings, even though the sorting (previous step) should be OK

bedtools bamtofastq -i ${name}.no_org.nodup.sorted.bam -fq ${name}_no_org_r1.fastq -fq2 ${name}_no_org_r2.fastq

rm ${name}.no_org.nodup.sorted.bam

# gzip (removes .fastq files; keeps .fastq.gz files)

gzip ${name}_no_org_r*.fastq

done

#
# run with
#
# ./script_name.sh 1> output_file.txt 2>/dev/null
#
#
#
###########################################



