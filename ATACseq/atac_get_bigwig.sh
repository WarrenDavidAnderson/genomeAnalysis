
# /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/preadip_bams/mm10_sob
# /scratch/wa3j/sob

# module information
module load genometools
module load samtools
module load ucsc-tools

# implement seqOutBias
gen=/nv/vol192/civeleklab/warren/MGlab/genomes/mm10/mm10.fa

for i in *atac.bam
do

bam=$(echo $i | awk -F"_sample_atac.bam" '{print $1}')

cat > sobScript_${bam}.sh <<EOF
#!/bin/bash
exec &> log_sob_${bam}.txt
echo ${bam}
nohup ./seqOutBias ${gen} ${i} --skip-bed --no-scale --bw=${bam}.bigWig --only-paired --shift-counts --read-size=38  &
EOF

# call a script for each core
echo sobScript_${bam}.sh
chmod 700 sobScript_${bam}.sh
nohup ./sobScript_${bam}.sh &

done


# just one
#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH -p largemem
#SBATCH --time=96:00:00
#SBATCH --mem=975G
#SBATCH --cpus-per-task=16
module load genometools
module load samtools
module load ucsc-tools
gen=/nv/vol192/civeleklab/warren/MGlab/genomes/mm10/mm10.fa
i=3T3_20min_rep1_sample_atac.bam
bam=$(echo $i | awk -F"_sample_atac.bam" '{print $1}')
./seqOutBias ${gen} ${i} --skip-bed --no-scale --bw=${bam}.bigWig --only-paired --shift-counts --read-size=38
