#!/bin/bash

# data directory
dir=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/makebedgraph
cd $dir

# path information
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/kentUtils/bin/linux.x86_64

# loop through each condition and merge bams
for ii in *rep1*scaled.bigWig
do

	# set the condition name, select files, start log
	name=$(echo $ii | awk -F"/" '{print $NF}' | awk -F"_scaled" '{print $1}')
	name=$(echo $name | awk -F"/" '{print $NF}' | awk -F"_rep" '{print $1}')		
	str=$(echo $ii | awk -F"/" '{print $NF}' | awk -F"_scaled" '{print $1}')
	str=$(echo $str | awk -F"/" '{print $NF}' | awk -F"rep1_" '{print $2}')	
	files=$(ls *${name}*${str}*scaled.bigWig*)
    exec &> logfile_${name}_${str}_conditionMerge.txt
    echo ${name}_${str}
	echo ${files}

	# merge the bigwigs
	bigWigMerge ${files} ${name}_${str}_scaled_merged.bedGraph

done
