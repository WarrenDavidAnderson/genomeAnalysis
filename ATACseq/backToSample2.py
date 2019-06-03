#############################

## code designed to remove reads that overlap in two input fastq files
## Warren Anderson, 8-2017

# input:
### fastq file with unique reads from sample species
### fastq file with unique reads from spike-in species

# outputs:
### fastq file with spike-in reads removed
### fastq file with overlap reads

# run instructions:
### python <name of this file>.py <input1> <input2> <output>
### example: python backToSample.py sample.fastq spike.fastq output.fastq
#############################

import sys
sys.path.insert(0, "/h4/t1/users/wa3j/software/py/biopython-1.70")
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
import sys


#############################
# specify input and output files
#############################

# specify name of input fastq files
input_ref = sys.argv[1] # sample species
input_spk = sys.argv[2] # spike-in species
fname_out = sys.argv[3] # output file name
#############################

# import read data
ref_reads = SeqIO.parse(open(input_ref, 'rU'), 'fastq') # SeqRecord object
spk_reads = SeqIO.parse(open(input_spk, 'rU'), 'fastq') # SeqRecord object

ref_id = [] 		# SEGUIDs
spk_id = [] 		# SEGUIDs
ref_out = [] 		# fastq reads
intersects = [] 		# fastq reads

# catalog spike-in reads (second input file)
for spk in spk_reads:
	spk_id.append( seguid(spk.seq) )

# catalog sample reads that do not overlap with spike-in reads
overlaps = 0
final = 0
spk = spk_id
output_handle1 = open(fname_out, "w")
output_handle2 = open("nomm_intersects.fastq", "w")
for ref in ref_reads: 
	cksum = seguid(ref.seq)
	if cksum not in spk:
		SeqIO.write(ref, output_handle1, 'fastq')
		final = final + 1
	if cksum in spk:
		SeqIO.write(ref, output_handle2, 'fastq')
		spk.remove( cksum ) 
		overlaps = overlaps + 1

output_handle1.close()
output_handle2.close()

# report overlaps
#print "\noverlaps detected:"
print "overlaps detected in no mismatch data" 
print overlaps
#print "\n"

# report final sample read count
print "\nread count after removing overlap reads" 
print final
#print "\n"


