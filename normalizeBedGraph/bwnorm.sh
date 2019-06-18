
# data directory
cd /media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/bigWig_20181001/browserBG

# path information
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/bedtools2/bin
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/genomeTools/bin
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/kentUtils/bin/linux.x86_64

# generate bedgraph files from un-normalized bigwigs
for ii in *.bigWig
do
cond=$(echo ${ii} | awk -F"_sample" '{print $1}')
strd=$(echo ${ii} | awk -F"_sample_" '{print $2}')
strd=$(echo ${strd} | awk -F".bigWig" '{print $1}')
bigWigToBedGraph ${ii} ${cond}_${strd}.bedgraph
done


############################################################
## R code - normalize each bedgraph
############################################################

library(dplyr)
library(DESeq2)
library(bigWig)

# import pro data
normfun = paste0("/run/user/1001/gvfs/smb-share:server=home1.virginia.edu,",
                 "share=cphgdesk/users/wa3j/MGlab/Analysis/Adipogenesis",
                 "/spike_norm/pro_normalization/pro_norm_functions.R")
source(normfun)

# TU annotation 
bed0 = read.table("TU_20181107.bed",stringsAsFactors=F,header=F)
names(bed0) = c("chr","start","end","gene","xy","strand")

# import pro data
bigWig.path = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                     "/July2018/bigWig_20181001/all_pro_sample_bigWig")
file.prefix = "3T3_"
file.suffix = "_sample_plus.bigWig"
min.length = 1
reads0 = get.counts.interval(bed=bed0, bigWig.path=bigWig.path,
                             file.prefix=file.prefix,
                             file.suffix=file.suffix,
                             min.length=min.length)

# estimate size factors
size_factors = estimateSizeFactorsForMatrix(reads0)

# normalize bg for each rep/strand
bgs = list.files(getwd(),".bedgraph")
for(ii in 1:length(bgs)){
cond = paste(strsplit(bgs[ii],"_")[[1]][2:3],collapse="_")
strd = strsplit(bgs[ii],"_")[[1]][4]
strd = strsplit(strd,".bed")[[1]][1]
dat0 = read.table(bgs[ii],header=F,stringsAsFactors=F)
sf = size_factors[names(size_factors)==cond]
dat0$V4 = dat0$V4 / sf
fname = paste0(cond,"_",strd,"_norm.bedGraph")
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
fpos=$(ls *.bedGraph | grep ${ii} | grep "plus")
fneg=$(ls *.bedGraph | grep ${ii} | grep "minus")
bedtools unionbedg -i ${fpos} |  awk -v OFS='\t' '{print $1, $2, $3, $4+$5+$6}' > pro_${ii}_plus.bedGraph
bedtools unionbedg -i ${fneg} |  awk -v OFS='\t' '{print $1, $2, $3, $4+$5+$6}' > pro_${ii}_minus.bedGraph
bedGraphToBigWig pro_${ii}_plus.bedGraph mm10.chrom.sizes pro_${ii}_plus.bigWig
bedGraphToBigWig pro_${ii}_minus.bedGraph mm10.chrom.sizes pro_${ii}_minus.bigWig
done




