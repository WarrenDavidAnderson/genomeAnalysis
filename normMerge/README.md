
# Scaling bigWigs/bedGraphs by DESeq2 size factors and merging

These analyses use Seqoutbias and UCSC tools. It is assumed that DESeq2 size factors have been identified. All analyses are completed in a strand-sensitive manner for PRO-seq data.

## bam_to_bigwig.sh

This script will loop through a directory of individual bam files generate bigWigs and bedGraphs, the scale the bedGraphs and convert to bigWigs for downstream analysis.

## merge_bigwigs.sh

This script will loop through conditions (e.g., "t0" and "20min") and merge the bigWigs for each condition into a resulting bedGraph file.

## bg_to_mergebigwig.sh

This script will go through each condition-merged bedGraph file and generate a cooresponding bigWig file.

## merge_all_bigwig.sh

This script will further merge all of the condition-merged files in a strand-specific manner. This was done to aggregate all pre-adipogenesis timeseries data.



