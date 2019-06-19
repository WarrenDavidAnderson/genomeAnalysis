
# data directory
cd /media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_bw_20181127/normdat

# path information
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/bedtools2/bin
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/genomeTools/bin
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/kentUtils/bin/linux.x86_64

# generate bedgraph files from un-normalized bigwigs
for ii in *.bigWig
do
cond=$(echo ${ii} | awk -F".bigWig" '{print $1}')
bigWigToBedGraph ${ii} ${cond}.bedgraph
done


############################################################
## R code - normalize each bedgraph
############################################################

library(dplyr)
library(DESeq2)
library(bigWig)

# import atac data
fun = paste0("/run/user/1001/gvfs/smb-share:server=home1.virginia.edu,",
			"share=cphgdesk/users/wa3j/MGlab/Analysis/Adipogenesis/",
			"ATAC_time_DEG/atac_norm_functions.R")
source(fun)

# atac peak annotation 
load("bed.map.RData")

# atac data mapped to peak regions
bigWig.path = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                     "/July2018/atac_bw_20181127/atac_all")
file.prefix = "3T3_"
file.suffix = ".bigWig"
min.length = 1
reads0 = get.counts.interval(bed=bed.map[,1:3], bigWig.path=bigWig.path,
                             file.prefix=file.prefix,
                             file.suffix=file.suffix,
                              min.length=min.length)

# estimate size factors
size_factors = estimateSizeFactorsForMatrix(reads0)

# normalize bg for each rep/strand
bgs = list.files(getwd(),".bedgraph")
for(ii in 1:length(bgs)){
cond = paste(strsplit(bgs[ii],"_")[[1]][2:3],collapse="_")
cond = strsplit(cond,".bedgraph")[[1]][1]
dat0 = read.table(bgs[ii],header=F,stringsAsFactors=F)
sf = size_factors[names(size_factors)==cond]
dat0$V4 = dat0$V4 / sf
fname = paste0(cond,"_norm.bedGraph")
write.table(dat0,fname,col.names=F,row.names=F,sep="\t",quote=F)
cat(paste0(ii / length(bgs)),"\n")
}


############################################################
## unix code to combine bedgraphs for each condition
############################################################

# get a list of conditions
cond=()
for ii in *.bedGraph
do
cc=$(echo ${ii} | awk -F"_rep" '{print $1}')
cond+=(${cc})
done

# merge bedgraphs for each unique condition, generage bigwigs as well
condu=($(printf "%s\n" "${cond[@]}" | sort -u | tr '\n' ' '))
for ii in "${condu[@]}"
do
echo ${ii}
file=$(ls *.bedGraph | grep ${ii})
bedtools unionbedg -i ${file} |  awk -v OFS='\t' '{print $1, $2, $3, $4+$5+$6}' > atac_${ii}.bedGraph
bedGraphToBigWig atac_${ii}.bedGraph mm10.chrom.sizes atac_${ii}.bigWig
done

############################################################
## add header to BGs for browser visualization
############################################################

# add header
for ii in *.bedGraph
do
mv ${ii} tmp
cond=$(echo ${ii} | awk -F"_" '{print $2}')
echo "track type=bedGraph name=${cond} color=0,0,0 altColor=0,0,0 alwaysZero=on visibility=full" > hdr
cat hdr tmp > ${ii}
rm hdr tmp
done

# compress
gzip *.bedGraph

