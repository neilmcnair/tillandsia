#!/bin/bash

# 5th script

# for comparing individual results with a 5 species comparative analysis

# Example command to run this script, from "Neil" folder on 8TB - "scripts/cluster_comp_by_read.sh | tee output.file"


set0='RepEx_results/5spp_paired_reads' # 5spp comparative run
set1='RepEx_results/BG9' # Tfasc
set2='RepEx_results/BG13' # Tleib
set3='RepEx_results/BG20' # Tprop
set4='RepEx_results/T_australis'
set5='RepEx_results/T_ionantha'


cutoff=0.1 # 0.25 seems like a reasonable value for 2spp, but 0.1 works better for 5spp

# i) It tries to find which clusters in a comparative run (RepEx) have lots of reads from clusters in individual runs


for comp in $(seq -f "%03g" 1 184) # need to check cluster count manually for each output dataset
do

# count the total number reads in the cluster in the comparative output

cltotreads=$(grep ">" ${set0}/seqclust/clustering/clusters/dir_CL0${comp}/reads.fasta | wc -l)

for first in $(seq -f "%03g" 1 246)
do

# count the reads in the comparative output which are present in the cluster from the first individual
# removing the first and last characters so that reads can be found from their pairs, where paired end data has been used

grep ">" ${set1}/seqclust/clustering/clusters/dir_CL0${first}/reads.fasta | sed 's/^.\(.*\).$/\1/' > temp1

readsfrom1stindiv=$(grep -f temp1 ${set0}/seqclust/clustering/clusters/dir_CL0${comp}/reads.fasta | wc -l)

# find the proportion of each comparative cluster that comes from each individual cluster

firstprop=$(echo "scale=2;$readsfrom1stindiv/$cltotreads" | bc)

# only print if at least a quarter of the comparative cluster corresponds to reads from individual clusters

if [ $(echo "scale=2;${firstprop} >= ${cutoff}" | bc) = 1 ]; then
  #echo "${comp} ${cltotreads} ${first} ${readsfrom1stindiv} ${second} ${readsfrom2ndindiv}"
  echo "${comp} "1st" ${first} ${firstprop}"
fi

# tidy up

rm temp1

done


for second in $(seq -f "%03g" 1 159) # manually add largest cluster number
do

# count the reads in the comparative output which are present in the cluster from the second individual

grep ">" ${set2}/seqclust/clustering/clusters/dir_CL0${second}/reads.fasta | sed 's/^.\(.*\).$/\1/' > temp2

readsfrom2ndindiv=$(grep -f temp2 ${set0}/seqclust/clustering/clusters/dir_CL0${comp}/reads.fasta | wc -l)

# find the proportion of each comparative cluster that comes from each individual cluster

secondprop=$(echo "scale=2;$readsfrom2ndindiv/$cltotreads" | bc)

# only print if at least a quarter of the comparative cluster corresponds to reads from individual clusters

if [ $(echo "scale=2;${secondprop} >= ${cutoff}" | bc) = 1 ]; then
  #echo "${comp} ${cltotreads} ${first} ${readsfrom1stindiv} ${second} ${readsfrom2ndindiv}"
  echo "${comp} "2nd" ${second} ${secondprop}"
fi

# tidy up

rm temp2

done

for third in $(seq -f "%03g" 1 202)
do

# count the reads in the comparative output which are present in the cluster from the second individual

grep ">" ${set3}/seqclust/clustering/clusters/dir_CL0${third}/reads.fasta | sed 's/^.\(.*\).$/\1/' > temp2

readsfrom2ndindiv=$(grep -f temp2 ${set0}/seqclust/clustering/clusters/dir_CL0${comp}/reads.fasta | wc -l)

# find the proportion of each comparative cluster that comes from each individual cluster

secondprop=$(echo "scale=2;$readsfrom2ndindiv/$cltotreads" | bc)

# only print if at least a quarter of the comparative cluster corresponds to reads from individual clusters

if [ $(echo "scale=2;${secondprop} >= ${cutoff}" | bc) = 1 ]; then
  #echo "${comp} ${cltotreads} ${first} ${readsfrom1stindiv} ${second} ${readsfrom2ndindiv}"
  echo "${comp} "3rd" ${third} ${secondprop}"
fi

# tidy up

rm temp2

done



for cl in $(seq -f "%03g" 1 199)
do

# count the reads in the comparative output which are present in the cluster from the second individual

grep ">" ${set4}/seqclust/clustering/clusters/dir_CL0${cl}/reads.fasta | sed 's/^.\(.*\).$/\1/' > temp # cheated here, as for ionantha need to say reads.fasta instead of .fas

readsfromindiv=$(grep -f temp ${set0}/seqclust/clustering/clusters/dir_CL0${comp}/reads.fasta | wc -l)

# find the proportion of each comparative cluster that comes from each individual cluster

prop=$(echo "scale=2;$readsfromindiv/$cltotreads" | bc)

# only print if at least a quarter of the comparative cluster corresponds to reads from individual clusters

if [ $(echo "scale=2;${prop} >= ${cutoff}" | bc) = 1 ]; then
  #echo "${comp} ${cltotreads} ${first} ${readsfrom1stindiv} ${second} ${readsfrom2ndindiv}"
  echo "${comp} "4th" ${cl} ${prop}"
fi

# tidy up

rm temp

done



for cl in $(seq -f "%03g" 1 202)
do

# count the reads in the comparative output which are present in the cluster from the fifth individual

grep ">" ${set5}/seqclust/clustering/clusters/dir_CL0${cl}/reads.fasta | sed 's/^.\(.*\).$/\1/' > temp # cheated here, as for juncea need to say reads.fasta instead of .fas

readsfromindiv=$(grep -f temp ${set0}/seqclust/clustering/clusters/dir_CL0${comp}/reads.fasta | wc -l) 

# find the proportion of each comparative cluster that comes from each individual cluster

prop=$(echo "scale=2;$readsfromindiv/$cltotreads" | bc)

# only print if at least the required proportion of the comparative cluster corresponds to reads from individual clusters

if [ $(echo "scale=2;${prop} >= ${cutoff}" | bc) = 1 ]; then
  #echo "${comp} ${cltotreads} ${first} ${readsfrom1stindiv} ${second} ${readsfrom2ndindiv}"
  echo "${comp} "5th" ${cl} ${prop}"
fi

# tidy up

rm temp

done

done


