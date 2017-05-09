
#!/usr/bin/env python

import sys, os
import argparse
from argparse import RawTextHelpFormatter
from collections import Counter #for counting characters in a string
from Bio import AlignIO




def main():



	alignment = AlignIO.read(args.input[0], "fasta")
	alignLen = alignment.get_alignment_length() #get length of alignment i.e. number of columns
	numTaxa = len(alignment) #count number of taxa

	if not args.output:
		output_file = "informativeOnly." + args.input[0]
		
	else:
	 	output_file = args.output
 


	print "\nProcessing an alignment with %s taxa and %s positions..." % (numTaxa, alignLen)



	listSeqs = {} #initialise list of sequences
	
	for k in alignment: 
		seq1 = str(k.seq).upper()
	
		if seq1 not in listSeqs:
			listSeqs[seq1] = k.description
		else:
			listSeqs[seq1] += "\t" + k.description
	
	noIdenticalOutfile = open("no_identicals." + output_file,"w")
	
	for k in listSeqs:
		noIdenticalOutfile.write(">" + listSeqs[k] + "\n" + k +"\n")
		
	print "\nIdentical sequences collapsed and written to %s" % (noIdenticalOutfile)
		
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
Remove identical sequences from a fasta alignment.
    
Input: Fasta
    
Requires: BLAST+, ACT and BioPython on your PATH
    ''', formatter_class=RawTextHelpFormatter)
        
    
    #takes in the input files
    parser.add_argument('input', nargs="+", action="store", help="Specify at least 2 input files")
    
    parser.add_argument("-o", "--output", action="store", default=False, help="Output file. Default ...")
    parser.add_argument("-c", "--collapse", action="store_true", default=False, help="Collapse identical sequences, after removing uninformative. A separate file is produced'")


    #unused arguments
    parser.add_argument("-f", "--format", action="store", help="Specify output format e.g. fasta or phylip'")

    args = parser.parse_args()


    main() #run main script

