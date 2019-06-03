
# Basic processing for PRO-seq fastq files

Run the following analysis files:

atacseq_alignment.sh  
atac_get_bigwig.sh  

## atacseq_alignment.sh

Note. This pipeline will isolate both sample and spike-in species reads through a process of iterative alignment and filtering.

-Align the data to the genome using bowtie2 (2.2.5, http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

-Filter and process the alignments using samtools (1.6, http://samtools.sourceforge.net/)

-Use bedtools to convert bam files to the fastq format for aligning to both sample and spike in species (2.19.1, https://bedtools.readthedocs.io/en/latest/)

-Custom python scripts are used for various tasks related to characterizing the alignment pipeline and isolating species-specific reads (pyShell.py, removeOverlap2.py, backToSample2.py)

## atac_get_bigwig.sh

Use seqOutBias to isolate bigWig files (v1.1.3, https://guertinlab.github.io/seqOutBias/)
