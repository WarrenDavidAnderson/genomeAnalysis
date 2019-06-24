#!/bin/bash

# data directory
dir=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/makebedgraph
cd $dir

# chrom sizes file
chromsizes=mm10.chrom.sizes

# path information
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/kentUtils/bin/linux.x86_64

# files
files=$(ls *scaled_merged*.bedGraph)

# loop through each condition and merge bams
for ii in $files
do

	# set the condition name, select files, start log		
	str=$(echo $ii | awk -F"/" '{print $NF}' | awk -F"_scaled" '{print $1}')
	str=$(echo $str | awk -F"/" '{print $NF}' | awk -F"_" '{print $2}')
	name=$(echo $ii | awk -F"/" '{print $NF}' | awk -F"_${str}" '{print $1}')

    exec &> logfile_${name}_${str}_conditionMergetoBigWig.txt
    echo ${name}_${str}


	# generate the bigwig
	bedGraphToBigWig ${ii} ${chromsizes} ${name}_${str}_scaled_merged.bigWig

done
