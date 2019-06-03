

cd /nv/vol192/civeleklab/warren/MGlab/PRO_WAFD/adipo_pro/bams_adip

module load samtools

# process plus and minus aligned reads separately
for i in *.bam*
do

name=$(echo $i | awk -F".bam" '{print $1}')

cat > proScript${name}.sh <<EOF
#!/bin/bash

samtools view -bh -F 20 ${name}.bam > ${name}_pro_plus.bam
samtools view -bh -f 0x10 ${name}.bam > ${name}_pro_minus.bam

EOF

# call a script for each core
echo calling proScript${name}.sh
chmod 700 proScript${name}.sh
nohup ./proScript${name}.sh &

done
