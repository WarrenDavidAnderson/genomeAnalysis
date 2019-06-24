#!/bin/bash

# get genome sizes
# wget https://raw.githubusercontent.com/igvteam/igv/master/genomes/sizes/mm10.chrom.sizes

# size factor file
sf_dir=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018
sf_file=${sf_dir}/bg_spikeNorm.txt

# data directories
bamdir=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/bams
datdir=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/makebedgraph

# key analysis files
reg=/media/wa3j/Seagate2/Documents/genomes/genes_mouse/Mus_musculus.GRCm38.87.body
gen=/media/wa3j/Seagate2/Documents/genomes/mm10/mm10.fa
chromsizes=mm10.chrom.sizes
pyscript=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/scripts/normalize_bedGraph.py

# path information
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/bowtie2-2.2.6
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/samtools-1.2
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/bedtools2/bin
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/genomeTools/bin
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/kentUtils/bin/linux.x86_64
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/seqOutBias_v1.1.1.bin.linux.x86_64

# set directory to bam file location
cd $bamdir

# specify bams of interest
bams=$(ls *.bam*)

# convert from bams to bedgraph/bigwig
for i in ${bams}
do

	# specify file name, strand, and set log
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_pro" '{print $1}')
    name=$(echo $name | awk -F"/" '{print $NF}' | awk -F"adi_" '{print $2}')
	str=$(echo $i | awk -F"/" '{print $NF}' | awk -F"pro_" '{print $2}')
    str=$(echo $str | awk -F"/" '{print $NF}' | awk -F".bam" '{print $1}')
    exec &> logfile_${name}_${str}.txt
    echo ${name}_${str}

	# identify the appropriate normalization factor
    fctr=$(grep ${name} ${sf_file} | awk -F" " '{print $2}')
    scale=$(echo ${fctr} |bc)
	scaletrue=$(bc <<< "scale=4 ; 1.0 / ${scale}")
    echo $scale
    echo $scaletrue

	# implement seqoutbias
    seqOutBias ${gen} ${i} --regions=${reg}.${str}.bed --skip-bed --no-scale --bw=${name}_${str}.bigWig --tail-edge --read-size=30
    echo seqOutbias_done

	# convert the bigwig to a bedgraph
    bigWigToBedGraph ${name}_${str}.bigWig ${name}_${str}.bg
    
	# normalize the bedgraph by the size factor
	echo normalizing
    python ${pyscript} -i ${name}_${str}.bg -s $scaletrue -o ${name}_${str}_scaled.bg

	# generate the normalized (individual) bedGraph file
    touch temp.txt
    echo "track type=bedGraph name=${name}_${str} color=0,0,255 altColor=0,0,255 alwaysZero=on visibility=full" >> temp.txt
    cat temp.txt ${name}_${str}_scaled.bg > ${name}_${str}_scaled.bedGraph

	# remove intermediate files    
	rm temp.txt
    rm ${name}_${str}.bg
    rm ${name}_${str}_scaled.bg
    rm *tmp

	# convert the scaled bedgraph to a bigwig
    bedGraphToBigWig ${name}_${str}_scaled.bedGraph ${chromsizes} ${name}_${str}_scaled.bigWig

	# move bigwig and bedgraph to the data directory
	mv *.bedGraph* ${datdir}
	mv *.bigWig* ${datdir}
    
done



