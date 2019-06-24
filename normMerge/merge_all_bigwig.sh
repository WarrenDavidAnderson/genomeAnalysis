

# data directory
dir=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/makebedgraph/all_preadip
cd $dir

# chrom sizes file
chromsizes=mm10.chrom.sizes

# path information
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/kentUtils/bin/linux.x86_64

# merge plus
plusfiles=$(ls *plus*)
bigWigMerge ${plusfiles} preadip_plus_scaled_merged.bedGraph
bedGraphToBigWig preadip_plus_scaled_merged.bedGraph ${chromsizes} preadip_plus_scaled_merged.bigWig

# merge minus
minusfiles=$(ls *minus*)
bigWigMerge ${minusfiles} preadip_minus_scaled_merged.bedGraph
bedGraphToBigWig preadip_minus_scaled_merged.bedGraph ${chromsizes} preadip_minus_scaled_merged.bigWig
