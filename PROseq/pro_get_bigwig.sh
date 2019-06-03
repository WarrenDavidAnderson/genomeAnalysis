

#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH -p standard
#SBATCH --time=48:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=4

# /nv/vol192/civeleklab/warren/MGlab/PRO_WAFD/adipo_pro/bams_adip

# path information
module load genometools
module load samtools
module load ucsc-tools

# implement seqOutBias
gen=/nv/vol192/civeleklab/warren/MGlab/genomes/mm10/mm10.fa

for i in *.bam
do

bam=$(echo $i | awk -F"_sample" '{print $1}')
strand=$(echo $i | awk -F"sample_pro_" '{print $2}')
strand=$(echo $strand | awk -F".bam" '{print $1}')

exec &> logfile_${bam}_${strand}.txt
./seqOutBias ${gen} ${i} --no-scale --bw=${bam}_${strand}.bigWig --tail-edge --read-size=30 --skip-bed 

done
