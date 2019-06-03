#############################
## code designed to identify, document, and remove putative PRC duplicates
## Warren Anderson, 7-2017

# input:
### fastq file obtained from sequencing with unique molecular identifier barcodes
##### note that multiple instances of identical sequences correspond to PCR duplicates
##### originating from the same original RNA

# outputs:
### fastq file for unique reads
### fastq file for non-unique (duplicated) reads
### txt file for counts of unique and non-unique read counts

# run instructions:
### python <name of this file>.py <name of input fastq file>
### example: python deduplicator.py sample.fastq
#############################

import sys
sys.path.insert(0, "/h4/t1/users/wa3j/software/py/biopython-1.70")
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid

#############################
# specify input file name
#############################

# specify name of an input fastq file
if (len(sys.argv) == 3):
	inputFile = sys.argv[1] 
elif (len(sys.argv) == 1):   
  print "Please supply a fastq filename and output (i.e. python script.py in.fastq out.fastq)"   
  sys.exit()
#############################



#############################
# specify output file names
#############################

# establish sample identifier 
# input sample name with ".fastq" removed
sample_id = inputFile.split('.fastq')[0]

# file for unique reads
deduplicated = sys.argv[2]
#############################


# import the SeqRecord indicator (Seq object with identifiers)
data0 = SeqIO.parse(inputFile, 'fastq')


##################################################################
### determine which reads are duplicated
### save duplicated and unique reads
### run a basic error check
##################################################################

##### function to find putative PCR duplicates
def duplicated_molecID(reads):
	
	### loop through each read in the SeqRecord (rec)
	compare = []
	output_handle = open(deduplicated, "w")
	for rec in reads:

		# save all reads and SEGUIDs for non-duplicated sequences
		if (seguid(rec.seq) != compare):
			SeqIO.write(rec, output_handle, "fastq")

		compare = seguid(rec.seq)

	### end loop

	output_handle.close()

###### end function
###### function to find putative PCR duplicates
##################################################################


##################################################################
### output files with reads
### conditioned on failure to detect code error
##################################################################

duplicated_molecID(data0)
##################################################################



