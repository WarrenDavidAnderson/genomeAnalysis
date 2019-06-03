

pth=/m/civeleklab/civeleklab/warren/MGlab/PRO_WAFD/adipo_pro
genome_sample=/m/civeleklab/civeleklab/warren/MGlab/genomes/mm10/mm10
genome_spike=/m/civeleklab/civeleklab/warren/MGlab/genomes/dm6/dm6
fastx=/h4/t1/users/wa3j/software/fastx
cutadapt=/h4/t1/apps/seqanal/cutadapt-1.8.1-py2.7/bin/cutadapt
PATH=$PATH:/h4/t1/apps/seqanal/cutadapt-1.8.1-py2.7/bin
PATH=$PATH:/h4/t1/users/wa3j/software/fastx
PATH=$PATH:/h4/t1/apps/seqanal/bowtie2-2.2.5
PATH=$PATH:/h4/t1/apps/seqanal/samtools-1.2
PATH=$PATH:/h4/t1/apps/seqanal/bedtools/bin
adapter_seq=TGGAATTCTCGGGTGCCAAGG
minlen=23
lread=9
rread=39

# set path to directory for the data
cd $pth

#####################################################################
## process each fastq file
## loop through multiple cores
#####################################################################

for i in *.fastq*
do

# set path to directory for the data
cd $pth

# basic file name, global variable
n0=$(echo $i | awk -F".fastq" '{print $1}')
n1=${n0}.fastq

# create a directory for this sample
# move the sample data to this directory for all subsequent processing
mkdir ${n0}
cp -t ${n0} ${n1} deduplicator.py removeOverlap2.py backToSample2.py pyShell.py
cd ${n0}

#####################################################################
## name key for intermediate files
####################################################################

# n1: input fastq base file name
# n2: input fastq file following removal of PCR duplicates, used for downstream processing
# n3: no longer used - NULL
# n4: trimmed fastq file
# n5: reverse-complemented fastq.gz, this will be aligned to both the sample and spike-in genomes
# n6: initial alignment to sample species	
# n7: sorted and MAPQ-filtered bam
# n8: sorted and MAPQ-filtered bam converted to a fastq (primary sample alignment)
# n8b: sorted and MAPQ-filtered sample fastq with intersect reads removed
# n9: unique sample reads aligned to spike-in species
# n10: sample reads aligned to spike-in species, sorted and MAPQ-filtered bam
# n11: sample reads aligned to spike-in species, sorted and MAPQ-filtered fastq (intersect)
# n12: overlap reads aligned to the sample with no mismatches
# n13: overlap reads aligned to the spike-in with no mismatches
# n14: MAPQ-filtered overlap reads mapped to the sample species with no mismatches - remove these reads from the spike-in
# n15: MAPQ-filtered overlap reads mapped to the spike-in species with no mismatches - remove these reads from the sample
# n16a: fastq of overlap reads mapped to the sample species with no mismatches
# n16b: fastq of overlap reads mapped to the spike-in species with no mismatches
# n17: fastq file of unique sample reads with spike-in reads removed - this will be re-aligned and used for down stream analysis
# n18: primary spike-in alignment of the PCR de-duplicated and fastx processed reads
# n19: primary spike-in alignment, MAPQ-filtered bam
# n20: fastq of uniquely mapped spike-in reads
# n20b: fastq of uniquely mapped spike-in reads with intersect removed
# n21: spike-in reads with sample reads removed
####################################################################

#####################################################################
## generate a script for each core
#####################################################################

cat > proScript${n0}.sh <<EOF
#!/bin/bash

##################################
## initial processing
## (1) clipping adapters
## (2) remove PCR duplicates
## (3) trimming sequences
## (4) reverse complement
##################################

# initialize log file
exec &> log_${n0}.txt
echo $n0
echo ""

# basic counts
echo "reads in original fastq"
echo $(python -c "print $(wc -l ${n1} | awk '{print $1}')/4.0" | bc) 
echo ""

# clipping the original fastq to quantify adapter-adapter ligation products
# print number of reads in clipped file to log
# note. fraction of non-adapter ligation product reads = clipped / total
${cutadapt} -m ${minlen} -a ${adapter_seq} -o ${n0}.clipped.fastq ${n1}
echo "total lines in clipped fastq (reads*4)"
python pyShell.py ${n0}.clipped.fastq
echo ""

# exclude molecular barcode repeats
# first sort the data
# n2: fastq file following removal of PCR duplicates
paste - - - - < ${n0}.clipped.fastq | sort --parallel 4 --buffer-size 5G --stable -t $'\t' -k2,2 | tr '\t' '\n' > sorted.fastq
python deduplicator.py sorted.fastq n2.fastq
echo "total lines in post-clipped deducplicated fastq (reads*4)"
python pyShell.py n2.fastq
echo ""
rm sorted.fastq

# fastx trimming
# n4: trimmed fastq
# note. fraction of reads present post-trimming = trimmed / total 
${fastx}/fastx_trimmer -Q 33 -f ${lread} -l ${rread} -i n2.fastq -o n4.fastq 
echo "total lines in trimmed fastq (reads*4)"
python pyShell.py n4.fastq
echo ""
rm ${n0}.clipped.fastq n2.fastq

# fastx reverse compliment
# n5: reverse-complemented fastq.gz
${fastx}/fastx_reverse_complement -Q 33 -z -i n4.fastq -o n5.fastq.gz 
rm n4.fastq

####################################################################
## isolate sample-specific reads for analysis
## (1) identify reads uniquely mapped to the sample genome and the spike-in genome 
## (2) re-align these overlapping to both genomes without mismatches to the seed
## (3) remove reads that aligned to the spike-in with no mismatches
## (4) re-align these reads to the sample genome for subsequent processing
####################################################################

# alignment - sample species
# n6: initial alignment to sample species
echo "initial alignment to sample species"
bowtie2 -p 4 -x ${genome_sample} -U n5.fastq.gz -S n6.sam
echo ""

# sort sam file, filter based on MAPQ=10, and generate bam file
# n7: sorted and MAPQ-filtered bam
samtools view -b -q 10 n6.sam | samtools sort - n7 
echo "number of sample MAPQ-filtered reads"
samtools view -c n7.bam
echo ""
rm n6.sam

# convert the bam file into a fastq file for re-alignment to the 'spike-in' species
# n8: sorted and MAPQ-filtered bam converted to a fastq - spike-in reads will be removed from this file
bedtools bamtofastq -i n7.bam -fq n8.fastq
rm n7.bam

# align processed read data to 'spike-in' species
# here we are aligning de-duplicated and 'uniquely' mapped reads from alignment to the sample species
# sort sam file and generate bam file
# n9: unique sample reads aligned to spike-in species
# n10: sample reads aligned to spike-in species, sorted and MAPQ-filtered bam
# n11: sample reads aligned to spike-in species, sorted and MAPQ-filtered fastq
echo "align sample-aligned and MAPQ-filtered data to 'spike-in' species"
bowtie2 -x ${genome_spike} -U n8.fastq -S n9.sam
echo ""	
samtools view -b -q 10 n9.sam | samtools sort - n10 
echo "number of sample MAPQ-filtered reads aligned to spike-in (intersect - doubly unique)"
samtools view -c n10.bam
echo ""
bedtools bamtofastq -i n10.bam -fq n11.fastq 
rm n9.sam n10.bam

# remove doubly unique (intersect) reads from sample
# n8b: sorted and MAPQ-filtered fastq with intersect reads removed
#inputFile1=n8.fastq  	# input 1 to python script - reads will be removed from this file
#inputFile2=n11.fastq	# input 2 to python script - these reads will be removed
#fileOut=n8b.fastq 	# output to python script
echo "remove intersect from sample"
python removeOverlap2.py n8.fastq n11.fastq n8b.fastq
echo ""
echo "line count after removing overlap reads"
python pyShell.py n8b.fastq
rm n8.fastq

# realign thes doubly unique intersect reads to the sample species 
# do not allow mismatches
# uniquely aligned reads should be removed from the spike-in reads
# n12: overlap reads aligned to the sample with no mismatches
echo ""
echo "align these doubly overlap to the sample species without mismatches" 
bowtie2 --no-1mm-upfront -x ${genome_sample} -U n11.fastq -S n12.sam
echo "" 

# realign these doubly unique to the spike-in species without mismatches
# uniquely aligned reads should be removed from the sample reads 
# n13: overlap reads aligned to the spike-in with no mismatches
echo "align these doubly overlap to the spike-in species without mismatches"
bowtie2 --no-1mm-upfront -x ${genome_spike} -U n11.fastq -S n13.sam
echo ""

# extract unique alignments in the sample and spike-in no mismatch files
# n14: MAPQ-filtered overlap reads mapped to the sample species with no mismatches
# n15: MAPQ-filtered overlap reads mapped to the spike-in species with no mismatches
samtools view -b -q 10 n12.sam | samtools sort - n14 # sample
samtools view -b -q 10 n13.sam | samtools sort - n15 # spike-in
echo "reads from no missmatch sample alignment of intersect"
samtools view -c n14.bam
echo ""
echo "reads from no missmatch spike-in alignment of intersect"
samtools view -c n15.bam
echo "" 
rm n12.sam n13.sam

# remove reads from the sample that align without mismatches to the spike-in
# n16a: fastq of overlap reads mapped to the sample species with no mismatches
# n16b: fastq of overlap reads mapped to the spike-in species with no mismatches
# n17: fastq file of unique sample reads with spike-in reads removed - this will be re-aligned and used for downstream analysis
# inputFile1=n11.fastq  	# input 1 to python script - reads will be removed from this file
# inputFile2=n16.fastq		# input 2 to python script - these reads will be removed 
# fileOut=backtosample 		# output to python script 
bedtools bamtofastq -i n14.bam -fq n16a.fastq
bedtools bamtofastq -i n15.bam -fq n16b.fastq
echo "remove spike in no-mm from sample nomm to get reads that shpuld go back to the sample"
python backToSample2.py n16a.fastq n16b.fastq backtosample.fastq
echo ""
cat n8b.fastq backtosample.fastq > n17.fastq
rm n8b.fastq backtosample.fastq n14.bam n15.bam

# re-align final sample reads and create a properly named bam file
echo "final alignment for sample data"
bowtie2 -x ${genome_sample} -U n17.fastq -S sample_output.sam
echo ""
samtools view -b -q 10 sample_output.sam | samtools sort - ${n0}_sample
echo "final sample read count"
samtools view -c ${n0}_sample.bam 
echo ""
rm n17.fastq sample_output.sam
####################################################################

####################################################################
## determine the count of spike-in-specific reads 
## (1) find reads that map uniquely to the spike-in genome
## (2) remove overlapping reads that aligned to the sample with no mismatches
## (3) count remaining reads
####################################################################

# alignment - spike in species
# n18: primary spike-in alignment of the PCR de-duplicated and fastx processed reads
# n19: primary spike-in alignment, MAPQ-filtered bam
echo "spike-in alignment"
bowtie2 -x ${genome_spike} -U n5.fastq.gz -S n18.sam
echo ""	
samtools view -b -q 10 n18.sam | samtools sort - n19
rm n18.sam n5.fastq.gz

echo "number of spike-in MAPQ-filtered reads"
samtools view -c n19.bam
echo ""

# remove doubly unique (intersect) reads from spike-in
# n20: fastq of uniquely mapped spike-in reads
# n20b: fastq of uniquely mapped spike-in reads with intersect removed
#inputFile1=n20.fastq  	# input 1 to python script - reads will be removed from this file
#inputFile2=n11.fastq	# input 2 to python script - these reads will be removed
#fileOut=n8_new.fastq 	# output to python script
bedtools bamtofastq -i n19.bam -fq n20.fastq # spike-in reads
echo "remove intersect from spike-in"
python removeOverlap2.py n20.fastq n11.fastq n20b.fastq
echo ""
rm n19.bam n20.fastq n11.fastq

# remove reads from the spike-in that align without mismatches to the sample
# n21: spike-in reads with sample reads removed
# inputFile1=n16b.fastq  	# input 1 to python script - reads will be removed from this file
# inputFile2=n16a.fastq		# input 2 to python script - these reads will be removed 
# fileOut=backtosample 		# output to python script
echo "remove sample no-mm reads from spike nomm to get reads that should go back to the spike" 
python backToSample2.py n16b.fastq n16a.fastq backtospike.fastq
echo ""
cat n20b.fastq backtospike.fastq > n21.fastq
rm *n16* backtospike.fastq n20b.fastq

# re-align final spike-in reads and create a properly named bam file
echo "final alignment to the spkike-in"
bowtie2 -x ${genome_spike} -U n21.fastq -S spike_output.sam
echo ""
samtools view -b -q 10 spike_output.sam | samtools sort - ${n0}_spikeIn
echo "final spike-in read count"
samtools view -c ${n0}_spikeIn.bam
rm n21.fastq spike_output.sam

EOF

# call a script for each core
echo calling proScript${n0}.sh
chmod 700 proScript${n0}.sh
nohup ./proScript${n0}.sh &

done

