#############################

## code designed to remove reads that overlap in two input fastq files
## Warren Anderson, 8-2017

# input:
### fastq file with unique reads from sample species
### fastq file with unique reads from spike-in species

# outputs:
### fastq file with spike-in reads removed

# run instructions:
### python <name of this file>.py <input1> <input2> <output>
### example: python removeOverlap.py sample.fastq spike.fastq output.fastq
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
if (len(sys.argv) == 4):
	input_ref = sys.argv[1] # sample species
        input_spk = sys.argv[2] # spike-in species
        fname_out = sys.argv[3] # output file name 
elif (len(sys.argv) < 4):   
  print "Please supply a input and output filenames (i.e. python script.py in1.fastq in2.fastq out.fastq)"   
  sys.exit()

#############################

# import read data
ref_reads = SeqIO.parse(open(input_ref, 'rU'), 'fastq') # SeqRecord object
spk_reads = SeqIO.parse(open(input_spk, 'rU'), 'fastq') # SeqRecord object

ref_id = [] 		# SEGUIDs
spk_id = [] 		# SEGUIDs
ref_out = [] 		# fastq reads

# catalog spike-in reads (second input file)
for spk in spk_reads:
	spk_id.append( seguid(spk.seq) )

# catalog sample reads that do not overlap with spike-in reads
overlaps = 0
spk = spk_id
output_handle = open(fname_out, "w")
for ref in ref_reads: # loop through each sample read
	cksum = seguid(ref.seq) # read ID
	if cksum not in spk: # check if the sample read is NOT in the set of spike reads
                SeqIO.write(ref, output_handle, "fastq")
	if cksum in spk: # check if the sample read is in the set of spike reads
                spk.remove( cksum )
		overlaps = overlaps + 1

output_handle.close()

# report overlaps
#print "\noverlaps detected:"
print "overlaps detected" 
print overlaps
#print "\n"










