#!/bin/sh

# 4th script

# Plan:

# i) from clustered_reads.fasta, take the read names we want

# ii) trim the suffix 1 or 2

# iii) deduplicate

# iv) from sequences.fasta, take the reads we want

names='BG9 BG13 BG20 T_australis T_ionantha'

for name in $names
do

# to put all clustered reads into one file -> cat ../RepEx_results_from_Cube/BG4_0.1x/seqclust/clustering/clusters/dir_CL*/reads.fas > BG4_clustered_reads.fasta

grep ">" $REPEX/RepEx_inputs/${name}_clustered_reads.fasta | sed 's/^......//' | sed 's/.\{1\}$//' | sort | uniq > badger

grep -A 1 -F -f badger $REPEX/RepEx_results_from_Cube/${name}_0.1x/seqclust/sequences/sequences.fasta > $REPEX/RepEx_inputs/${name}_repaired_clust_reads.fasta

rm badger

done

# still need to manually reinsert 'species_name' or whatever into the read name - use Geany
# also use Geany to remove double dashes where they occur

# need to get missing "half-pairs from /seqclust/reads/reads.fasta*"

###

grep ">" BG9_clustered_reads.fasta | sed 's/^.//' | sed 's/.\{1\}$//' | sort | uniq > badger

grep -A 1 -F -f badger BG9_input_reads.fasta > BG9_repaired_clust_reads.fasta

#---------------------------------------------------------------------------

# to create a single file with all input reads - FALSE

cat seqclust/reads/reads.fasta* > species_input_reads.fasta

# only need to

cp seqclust/reads/reads.fasta > species_input_reads.fasta

# to deduplicate the input reads (need to change line length from default 60) - no longer needed

seqkit rmdup -n -w 120 < species_input_reads.fasta > species_dedup_input_reads.fasta

# to extract readnames (directionless) from clustered reads

grep ">" BG9_clustered_reads.fasta | sed 's/^.//' | sed 's/.\{1\}$//' | sort | uniq > badger

# to re-create ("re-pair") clustered reads file, but with missing half-pairs restored

grep -A 1 -F -f badger BG9_input_reads.fasta > BG9_repaired_clust_reads.fasta

###



