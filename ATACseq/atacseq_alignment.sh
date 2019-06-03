
# basic information
cell=3T3
pth=/m/civeleklab/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3
sampleGenome=/m/civeleklab/civeleklab/warren/MGlab/genomes/mm10/mm10
spikeGenome=/m/civeleklab/civeleklab/warren/MGlab/genomes/dm6/dm6
cd ${pth}

# key directories and path information
PATH=$PATH:/h4/t1/users/wa3j/software/fastx
PATH=$PATH:/h4/t1/apps/seqanal/bowtie2-2.2.5
PATH=$PATH:/h4/t1/apps/seqanal/samtools-1.6
PATH=$PATH:/h4/t1/apps/seqanal/bedtools/bin

# loop through each file for a given cell-type designation
for fq in *${cell}*PE1.fastq.gz
do

cd ${pth}

# annotate the sample name, create a new directory, and move the data
name=$(echo ${fq} | awk -F"_PE1.fastq.gz" '{print $1}')
mkdir ${name}
cp *${name}*.gz ${name}
cp -t ${name} backToSample2.py removeOverlap2.py pyShell.py
cd ${name}


####################################################################
# write a shell script to process each sample in the loop separately
cat > atacScript_${name}.sh <<EOF
#!/bin/bash

# initialize log file
exec &> log_${name}.txt
echo ${name}
echo ""

# basic counts
printf "reads in original fastq \n"
echo $(python -c "print $(zcat ${name}_PE1.fastq.gz | wc -l | awk '{print $1}')/4.0" | bc) 
echo $(python -c "print $(zcat ${name}_PE2.fastq.gz | wc -l | awk '{print $1}')/4.0" | bc) 


##########################################
## isolate sample-specific reads
##########################################

# align to the sample genome
printf "\ninitial alignment to sample species \n"
bowtie2 --maxins 140 -x ${sampleGenome} \
-1 ${name}_PE1.fastq.gz -2 ${name}_PE2.fastq.gz -S ${name}_smp.sam

# sort the sam file, filter based on MAPQ=10, remove duplicates, and generate a bam file
samtools view -b -q 10 ${name}_smp.sam | samtools sort -n - | samtools fixmate -m - - | \
samtools sort - | samtools markdup -r - ${name}_smp_rmdup.bam
printf "\nnumber of sample MAPQ-filtered and deduplicated reads \n"
samtools view -c ${name}_smp_rmdup.bam
rm ${name}_smp.sam

# convert the bam file into a fastq file for re-alignment to the 'spike-in' species
# spike-in reads will be removed from the generated fastq files
samtools sort -n ${name}_smp_rmdup.bam -o ${name}_smp_sort.bam
bedtools bamtofastq -i ${name}_smp_sort.bam -fq ${name}_smp_r1.fastq -fq2 ${name}_smp_r2.fastq
rm ${name}_smp_rmdup.bam ${name}_smp_sort.bam

# align processed read data to 'spike-in' species
# here we are aligning de-duplicated and 'uniquely' mapped reads from alignment to the sample species
# the reads aligned here are the 'intersect' reads  
# these reads will be removed from the sample and re-aligned with no mismatches to seed
printf "\nalign sample-aligned and MAPQ-filtered data to 'spike-in' species \n"
bowtie2 --maxins 140 -x ${spikeGenome} \
-1 ${name}_smp_r1.fastq -2 ${name}_smp_r2.fastq -S ${name}_spk.sam
samtools view -b -q 10 ${name}_spk.sam | samtools sort -n - | samtools fixmate -m - - | \
samtools sort - | samtools markdup -r - ${name}_spk_rmdup.bam
printf "\nnumber of sample MAPQ-filtered and deduplicated reads \n"
samtools view -c ${name}_spk_rmdup.bam
samtools sort -n ${name}_spk_rmdup.bam -o ${name}_spk_sort.bam
bedtools bamtofastq -i ${name}_spk_sort.bam -fq ${name}_spk_r1.fastq -fq2 ${name}_spk_r2.fastq
rm ${name}_spk.sam ${name}_spk_rmdup.bam ${name}_spk_sort.bam

# remove doubly unique (intersect) reads from sample
# input 1 to python script - reads will be removed from this file
# input 2 to python script - these reads will be removed
# input 3 to python script - output file
printf "\nremove intersect from sample \n"
python removeOverlap2.py ${name}_smp_r1.fastq ${name}_spk_r1.fastq ${name}_smp_r1_spkrem.fastq
python removeOverlap2.py ${name}_smp_r2.fastq ${name}_spk_r2.fastq ${name}_smp_r2_spkrem.fastq
printf "\nline count after removing overlap reads \n"
python pyShell.py ${name}_smp_r1_spkrem.fastq
python pyShell.py ${name}_smp_r2_spkrem.fastq
rm ${name}_smp_r1.fastq ${name}_smp_r2.fastq 

# realign thes doubly unique intersect reads to the sample species 
# do not allow mismatches
printf "\nalign these doubly overlap to the sample species without mismatches \n" 
bowtie2 --maxins 140 --no-1mm-upfront -x ${sampleGenome} \
-1 ${name}_spk_r1.fastq -2 ${name}_spk_r2.fastq -S ${name}_smp_spkrem_smp.sam

# realign these doubly unique intersect reads to the spike-in species without mismatches
# uniquely aligned reads should be removed from the sample reads 
printf "\nalign these doubly overlap to the spike-in species without mismatches \n"
bowtie2 --maxins 140 --no-1mm-upfront -x ${spikeGenome} \
-1 ${name}_spk_r1.fastq -2 ${name}_spk_r2.fastq -S ${name}_smp_spkrem_spk.sam

# extract unique alignments in the sample no mismatch file
samtools view -b -q 10 ${name}_smp_spkrem_smp.sam | samtools sort -n - | \
samtools fixmate -m - - | samtools sort - | \
samtools markdup -r - ${name}_smp_spkrem_smp_rmdup.bam
printf "\nnumber of sample no-mm reads \n"
samtools view -c ${name}_smp_spkrem_smp_rmdup.bam
samtools sort -n ${name}_smp_spkrem_smp_rmdup.bam -o ${name}_smp_spkrem_smp_sort.bam
bedtools bamtofastq -i ${name}_smp_spkrem_smp_sort.bam \
-fq ${name}_smp_spkrem_smp_r1.fastq -fq2 ${name}_smp_spkrem_smp_r2.fastq
rm ${name}_smp_spkrem_smp.sam ${name}_smp_spkrem_smp_rmdup.bam
rm ${name}_smp_spkrem_smp_sort.bam
rm ${name}_spk_r1.fastq ${name}_spk_r2.fastq

# extract unique alignments in the spike-in no mismatch file
samtools view -b -q 10 ${name}_smp_spkrem_spk.sam | samtools sort -n - | \
samtools fixmate -m - - | samtools sort - | \
samtools markdup -r - ${name}_smp_spkrem_spk_rmdup.bam
printf "\nnumber of spike-in no-mm reads \n"
samtools view -c ${name}_smp_spkrem_spk_rmdup.bam
samtools sort -n ${name}_smp_spkrem_spk_rmdup.bam -o ${name}_smp_spkrem_spk_sort.bam
bedtools bamtofastq -i ${name}_smp_spkrem_spk_sort.bam \
-fq ${name}_smp_spkrem_spk_r1.fastq -fq2 ${name}_smp_spkrem_spk_r2.fastq
rm ${name}_smp_spkrem_spk.sam ${name}_smp_spkrem_spk_rmdup.bam
rm ${name}_smp_spkrem_spk_sort.bam

# remove reads from the no-mm sample that align without mismatches to the spike-in
# input 1 to python script - reads will be removed from this file
# input 2 to python script - these reads will be removed 
# input 3 to python script - output file
printf "\nremove spike in no-mm from sample no-mm reads to get reads that should go back to the sample \n"
python backToSample2.py ${name}_smp_spkrem_smp_r1.fastq ${name}_smp_spkrem_spk_r1.fastq ${name}_smp_spkrem2_r1.fastq
python backToSample2.py ${name}_smp_spkrem_smp_r2.fastq ${name}_smp_spkrem_spk_r2.fastq ${name}_smp_spkrem2_r2.fastq

# add the sample-specific no-mm reads back to the sample
cat ${name}_smp_r1_spkrem.fastq ${name}_smp_spkrem2_r1.fastq > ${name}_smp_r1.fastq
cat ${name}_smp_r2_spkrem.fastq ${name}_smp_spkrem2_r2.fastq > ${name}_smp_r2.fastq
rm ${name}_smp_r1_spkrem.fastq ${name}_smp_spkrem2_r1.fastq
rm ${name}_smp_r2_spkrem.fastq ${name}_smp_spkrem2_r2.fastq

# re-align final sample reads and generate the output bam file
printf "\nfinal alignment for sample data \n"
bowtie2 --maxins 140 -x ${sampleGenome} -1 ${name}_smp_r1.fastq -2 ${name}_smp_r2.fastq -S ${name}_sample1.sam
samtools view -b -q 10 ${name}_sample1.sam | samtools sort -n - | \
samtools fixmate -m - - | samtools sort - | \
samtools markdup -r - ${name}_sample_atac.bam
printf "\nfinal number of sample reads \n"
samtools view -c ${name}_sample_atac.bam
rm ${name}_sample1.sam 


####################################################################
## isolate spike-in-specific reads 
####################################################################

# align to the spike genome
printf "\ninitial alignment to spike-in species \n"
bowtie2 --maxins 140 -x ${spikeGenome} \
-1 ${name}_PE1.fastq.gz -2 ${name}_PE2.fastq.gz -S ${name}_spk.sam

# sort the sam file, filter based on MAPQ=10, remove duplicates, and generate a bam file
samtools view -b -q 10 ${name}_spk.sam | samtools sort -n - | samtools fixmate -m - - | \
samtools sort - | samtools markdup -r - ${name}_spk_rmdup.bam
printf "\nnumber of spike-in MAPQ-filtered and deduplicated reads \n"
samtools view -c ${name}_spk_rmdup.bam
rm ${name}_spk.sam

# convert the bam file into a fastq file for further processing
samtools sort -n ${name}_spk_rmdup.bam -o ${name}_spk_sort.bam
bedtools bamtofastq -i ${name}_spk_sort.bam -fq ${name}_spk_r1.fastq -fq2 ${name}_spk_r2.fastq
rm ${name}_spk_rmdup.bam ${name}_spk_sort.bam

# align processed read data to 'sample' species
# here we are aligning de-duplicated and 'uniquely' mapped reads from alignment to the spike-in species
# the reads aligned here are the 'intersect' reads  
# these reads will be removed from the spike-in and re-aligned with no mismatches to seed
printf "\nalign spike-aligned and MAPQ-filtered data to 'sample' species \n"
bowtie2 --maxins 140 -x ${sampleGenome} \
-1 ${name}_spk_r1.fastq -2 ${name}_spk_r2.fastq -S ${name}_spk.sam
samtools view -b -q 10 ${name}_spk.sam | samtools sort -n - | samtools fixmate -m - - | \
samtools sort - | samtools markdup -r - ${name}_spk_rmdup.bam
printf "\nnumber of spike MAPQ-filtered and deduplicated reads \n"
samtools view -c ${name}_spk_rmdup.bam
samtools sort -n ${name}_spk_rmdup.bam -o ${name}_spk_sort.bam
bedtools bamtofastq -i ${name}_spk_sort.bam -fq ${name}_spkrm_r1.fastq -fq2 ${name}_spkrm_r2.fastq
rm ${name}_spk.sam ${name}_spk_rmdup.bam ${name}_spk_sort.bam

# remove doubly unique (intersect) reads from spike-in
# input 1 to python script - reads will be removed from this file
# input 2 to python script - these reads will be removed
# input 3 to python script - output file
printf "\nremove intersect from sample \n"
python removeOverlap2.py ${name}_spk_r1.fastq ${name}_spkrm_r1.fastq ${name}_spk_r1_smprem.fastq
python removeOverlap2.py ${name}_spk_r2.fastq ${name}_spkrm_r2.fastq ${name}_spk_r2_smprem.fastq
printf "\nline count after removing overlap reads \n"
python pyShell.py ${name}_spk_r1_smprem.fastq
python pyShell.py ${name}_spk_r2_smprem.fastq
rm ${name}_spkrm_r1.fastq ${name}_spkrm_r2.fastq
rm ${name}_spk_r1.fastq ${name}_spk_r2.fastq

# remove reads from the no-mm spike-in that align without mismatches to the sample
# input 1 to python script - reads will be removed from this file
# input 2 to python script - these reads will be removed 
# input 3 to python script - output file
printf "\nremove sample in no-mm from spike no-mm reads to get reads that should go back to the spike-in \n"
python backToSample2.py ${name}_smp_spkrem_spk_r1.fastq ${name}_smp_spkrem_smp_r1.fastq backtospike_r1.fastq
python backToSample2.py ${name}_smp_spkrem_spk_r2.fastq ${name}_smp_spkrem_smp_r2.fastq backtospike_r2.fastq
cat ${name}_spk_r1_smprem.fastq backtospike_r1.fastq > finalspike_r1.fastq
cat ${name}_spk_r2_smprem.fastq backtospike_r2.fastq > finalspike_r2.fastq
rm ${name}_smp_spkrem_spk_r1.fastq ${name}_smp_spkrem_smp_r1.fastq
rm ${name}_smp_spkrem_spk_r2.fastq ${name}_smp_spkrem_smp_r2.fastq

# re-align final spike-in reads and generate the output bam file
printf "\nfinal alignment for spike-in data \n"
bowtie2 --maxins 140 -x ${spikeGenome} -1 finalspike_r1.fastq -2 finalspike_r2.fastq -S ${name}_spike1.sam
samtools view -b -q 10 ${name}_spike1.sam | samtools sort -n - | \
samtools fixmate -m - - | samtools sort - | \
samtools markdup -r - ${name}_spike_atac.bam
printf "\nfinal number of spike reads \n"
samtools view -c ${name}_spike_atac.bam
rm *_r1.fastq *_r2.fastq *_smprem.fastq ${name}_spike1.sam nomm_intersects.fastq


EOF
# end of shell script to process each sample in the loop separately
####################################################################

# call a script for each core
echo calling atacScript_${name}.sh
chmod 700 atacScript_${name}.sh
nohup ./atacScript_${name}.sh &

done

