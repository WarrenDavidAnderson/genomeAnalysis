

################################################################
## get data files, perform basic processing
## just run this in R
################################################################

lib.loc = "/h4/t1/users/wa3j/software/R_libs"
.libPaths(lib.loc)
library(dplyr, lib.loc=lib.loc)
library(bigWig, lib.loc=lib.loc)
library(groHMM, lib.loc=lib.loc)
library(GenomicFeatures, lib.loc=lib.loc)

# get merged bam file
bam.file="preadip_merged.bam"
data <- readGAlignments(bam.file)
data_gr <- granges(data)
data_gr <- sort(data_gr)
expr <- keepStandardChromosomes(data_gr, pruning.mode="coarse")

# get gene annotations based on TSS/TSS identification
gene.ann0 = read.table("inferredcoords.bed",sep="\t",header=F,stringsAsFactors=F)
names(gene.ann0) = c("chr","start","end","gene","xy","strand")
gene.ann = makeGRangesFromDataFrame(gene.ann0[,-5], seqnames.field="chr",
                                    start.field="start",end.field="end",
                                    strand.field="strand",
                                    starts.in.df.are.0based=TRUE,
                                    keep.extra.columns=TRUE)
gene.ann <- sort(gene.ann)
gene.ann <- keepStandardChromosomes(gene.ann, pruning.mode="coarse")

# window parameters
Fp <- windowAnalysis(expr, strand="+", windowSize=50)
Fm <- windowAnalysis(expr, strand="-", windowSize=50)

save.image("groHMMdata.RData")


################################################################
## HMM parameter testing
## hmm.param.test.R
################################################################

# run code 
# nohup Rscript hmm.param.test.R &

lib.loc = "/h4/t1/users/wa3j/software/R_libs"
.libPaths(lib.loc)
library(dplyr, lib.loc=lib.loc)
library(bigWig, lib.loc=lib.loc)
library(groHMM, lib.loc=lib.loc)
library(GenomicFeatures, lib.loc=lib.loc)

# load initial data
load("groHMMdata.RData")

# set the number of cores for grohmm
options(mc.cores=getCores(60))

# specify parameters for variation
vars = c(10,15,20,25,30,40,50,60,80,100,150,200,300)
ltpr = c(-50,-100,-200,-300,-400,-500)
LtProbB = sapply(ltpr,function(x){rep(x,length(vars))}) %>% as.vector
UTS = rep(vars,length(ltpr))
tune <- data.frame(LtProbB=LtProbB, UTS=UTS)

# testing parameter sets
hmm.vars <- mclapply(seq_len(nrow(tune)), function(x) {
  hmm <- detectTranscripts(Fp=Fp, Fm=Fm, LtProbB=tune$LtProbB[x], 
                           UTS=tune$UTS[x], threshold=1)
  e <- evaluateHMMInAnnotations(hmm$transcripts, gene.ann)
  return(list(hmm=hmm, eval=e$eval))
}, mc.cores=getOption("mc.cores"),  mc.silent=FALSE)
names(hmm.vars) = apply(tune,1,function(x){paste0(x[2],"_",x[1])})

save(hmm.vars, file="grohmm_vars.RData")


################################################################
## mouse genome annotation
################################################################

# set mm10 genome bed
# wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
gnme0 = read.table("mm10.chrom.sizes",header=F,stringsAsFactors=F)
names(gnme0) = c("chr","len")
rem = c(grep("random",gnme0$chr), grep("chrUn",gnme0$chr))
mm10.bed = gnme0[-rem,]
names(mm10.bed) = c("chr","end")
write.table(mm10.bed,"mm10.bed",sep="\t",quote=F,col.names=F,row.names=F)

# sort the file
command2=paste('sort -k1,1 -k2,2n', 'mm10.bed', '> mm10.sorted.bed')
system(command2)


################################################################
## get sensitivity information
## hmm.param.sens.R
################################################################

# run code 
# nohup Rscript hmm.param.sens.R &

# load libraries
lib.loc = "/h4/t1/users/wa3j/software/R_libs"
.libPaths(lib.loc)
library(bigWig, lib.loc=lib.loc)
library(groHMM, lib.loc=lib.loc)
library(GenomicFeatures, lib.loc=lib.loc)
library(dplyr, lib.loc=lib.loc)

# load hmm data
load("grohmm_vars.RData")

# load bigwigs
bw.plus = load.bigWig("preadip_plus_merged.bigWig")
bw.minus = load.bigWig("preadip_minus_merged.bigWig")

# bed dir
bed.bin = "/h4/t1/apps/seqanal/bedtools/bin/"

# loop through hmms, check reads in complement, document metrics

eval.dat = c()

for(ii in 1:length(hmm.vars)){
  
  # get complement to hmm regions
  hmm.test = hmm.vars[[ii]]$hmm
  transcripts = hmm.test[["transcripts"]]
  if(is.null(transcripts)==TRUE){next}
  hmm.bed0 <- data.frame(chr=seqnames(transcripts),
                         start=start(transcripts)-1,
                         end=end(transcripts),
                         names=c(rep(".", length(transcripts))),
                         scores=c(rep(".", length(transcripts))),
                         strand=strand(transcripts))
  hmmp = hmm.bed0 %>% filter(strand=="+")
  hmmm = hmm.bed0 %>% filter(strand=="-")
  write.table(hmmp,"hmmp.bed",sep="\t",quote=F,col.names=F,row.names=F)
  write.table(hmmm,"hmmm.bed",sep="\t",quote=F,col.names=F,row.names=F)
  command1=paste('sort -k1,1 -k2,2n', 'hmmp.bed', '> hmmp.sorted.bed')
  command2=paste('sort -k1,1 -k2,2n', 'hmmm.bed', '> hmmm.sorted.bed')
  system(command1); system(command2)
  comm1 = paste0(bed.bin,"complementBed -i hmmp.sorted.bed -g mm10.sorted.bed > compp.bed")
  comm2 = paste0(bed.bin,"complementBed -i hmmm.sorted.bed -g mm10.sorted.bed > compm.bed")
  system(comm1); system(comm2)
  complplus = read.table("compp.bed",stringsAsFactors=F)
  complminus = read.table("compm.bed",stringsAsFactors=F)
  system(paste0("rm hmmp.bed hmmm.bed hmmp.sorted.bed hmmm.sorted.bed compp.bed compm.bed"))
  
  # map reads to the complement regions
  plus.bed = complplus %>% mutate(gene=".",xy=".",strand="+")
  minus.bed = complplus %>% mutate(gene=".",xy=".",strand="-")
  names(plus.bed)[1:3] = names(minus.bed)[1:3] = c("chr","start","end")
  cnt.plus = bed.region.bpQuery.bigWig(bw.plus, plus.bed)
  cnt.minus = bed.region.bpQuery.bigWig(bw.minus, minus.bed)
  
  # get count density metrics
  len.p = plus.bed$end - plus.bed$start
  len.m = minus.bed$end - minus.bed$start
  density.plus = cnt.plus / len.p
  density.minus = cnt.minus / len.m
  counts = c(cnt.plus, cnt.minus)
  densities = c(density.plus, density.minus)
  max_cnt = max(counts);    max_den = max(densities)
  mean_cnt = mean(counts);  mean_den = mean(densities)
  med_cnt = median(counts); med_den = median(densities)
  sum_cnt = sum(counts);    sum_den = sum(densities)
  
  # combine evaluation metrics
  pars0 = strsplit(names(hmm.vars)[ii],"_")[[1]]
  LtProbB = pars0[2]
  UTS = pars0[1]
  hmm.eval = hmm.vars[[ii]]$eval
  new = data.frame(LtProbB, UTS, hmm.eval, 
                   max_cnt, mean_cnt, med_cnt, sum_cnt,
                   max_den, mean_den, med_den, sum_den)
  eval.dat = rbind(eval.dat, new)
  
} # ii, hmm loop

save(eval.dat, file="hmmEval.RData")


################################################################
## process error information to select the best hmm
################################################################

lib.loc = "/h4/t1/users/wa3j/software/R_libs"
.libPaths(lib.loc)
library(dplyr, lib.loc=lib.loc)
library(groHMM, lib.loc=lib.loc)
library(GenomicFeatures, lib.loc=lib.loc)

load("hmmEval.RData")
load("grohmm_vars.RData") 

# sort data based on merge errors, dissiciation errors, and sums of reads outside of TUs
eval.dat.sorted1 = eval.dat[with(eval.dat, order(merged, dissociated, sum_cnt)),] %>% 
  dplyr::select(merged, dissociated, sum_cnt)
eval.dat.sorted2 = eval.dat[with(eval.dat, order(dissociated, merged, sum_cnt)),] %>% 
  dplyr::select(merged, dissociated, sum_cnt)
eval.dat.sorted3 = eval.dat[with(eval.dat, order(sum_cnt, dissociated, merged)),] %>% 
  dplyr::select(merged, dissociated, sum_cnt)
head(eval.dat.sorted1)
head(eval.dat.sorted2)
head(eval.dat.sorted3)

# look at quartiles for metrics of interest
quantile(eval.dat$merged)
quantile(eval.dat$dissociated)
quantile(eval.dat$sum_cnt)

# select the best hmm: third lowest merge error (2106, lowest quartile)
# lowest sum count (7968454), lowest quartile for dissociation error (114)
eval.dat[68,]
hmm.best = hmm.vars[["20_-500"]]$hmm

# convert to bed
transcripts = hmm.best[["transcripts"]]
hmm.bed0 <- data.frame(chr=seqnames(transcripts),
                       start=start(transcripts)-1,
                       end=end(transcripts),
                       names=c(rep(".", length(transcripts))),
                       scores=c(rep(".", length(transcripts))),
                       strand=strand(transcripts))

save(hmm.bed0, file="bestHMM.bed")

