
# Basic processing for PRO-seq fastq files

Run the following analysis files:

proseq_alignment.sh
pro_get_plusminus.sh
pro_get_bigwig.sh

## proseq_alignment.sh

Note. This pipeline will isolate both sample and spike-in species reads through a process of iterative alignment and filtering.

-Clip adapter sequences with cutadapt (v1.8.1, https://cutadapt.readthedocs.io/en/stable/)

-Exclude molecular barcode repeats (deduplicator.py)

-Trim the reads using the fastx trimmer (http://hannonlab.cshl.edu/fastx_toolkit/)

-Reverse complement the reads using fastx tools (http://hannonlab.cshl.edu/fastx_toolkit/)

-Align the data to the genome using bowtie2 (v2.2.5, http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

-Filter and process the alignments using samtools (v1.2, http://samtools.sourceforge.net/)

-Use bedtools to convert bam files to the fastq format for aligning to both sample and spike in species (2.19.1, https://bedtools.readthedocs.io/en/latest/)

-Custom python scripts are used for various tasks related to characterizing the alignment pipeline and isolating species-specific reads (pyShell.py, removeOverlap2.py, backToSample2.py)

## pro_get_plusminus.sh

Use samtools to separate plus and minus strand reads (v1.9, http://samtools.sourceforge.net/)

## pro_get_bigwig.sh

Use seqOutBias to isolate strand-specific bigWig files (v1.1.3, https://guertinlab.github.io/seqOutBias/)


