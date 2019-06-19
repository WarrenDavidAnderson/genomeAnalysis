
# Normalize individual bigWig/bedGraph files and merge replicates

Analysis files:

bwnorm_pro.sh  
bwnorm_atac.sh  
pro_norm_functions.R  
atac_norm_functions.R  

## bwnorm_pro.sh and bwnorm_atac.sh

The following analyses are implemented using DESeq2 (R), Bedtools, and UCSCtools. Note that all analyses are completed in a strand-specific manner for use with PRO-seq data (bwnorm_pro.sh) but not for ATAC-seq data (bwnorm_atac.sh).

- Generate a bedGraph file for each bigWig replicate file (bigWigToBedGraph)

- Load replicate bigWig files in R and identify DESeq2 size factors, based on an existing bed coordinate gene/peak annotation (e.g., TU_20181107.bed)

- Load the individual replicate bedGraphs into R and normalize by the identified size factors

- For each condition, aggregate bedGraphs (unionbedg) and generate corresponding bigWig files (bedGraphToBigWig)

## pro_norm_functions.R

This file includes an R function for combining raw reads from bigWig files into a matrix.

## atac_norm_functions.R

This file includes an R function for combining raw reads from bigWig files into a matrix.


