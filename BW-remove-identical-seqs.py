#!/usr/bin/env python

#20170915 - BW edit. Added the ability to compare the percent identity between 2 strings

import sys, os
import argparse
from argparse import RawTextHelpFormatter
from collections import Counter #for counting characters in a string
from Bio import AlignIO
from difflib import SequenceMatcher


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

		#loop through listSeqs, if higher percent identity to a seq already in listSeqs. Change "aboveThreshold" to False
		aboveThreshold = False
		match = ""
		for seq2 in listSeqs:
			if similar(seq1, seq2) > args.threshold:
					aboveThreshold = True
					match = seq2 #store the sequence that matched. This will be the representative sequence.
			else:
				next

		#If there is a sequence > threshold, append sequence ID to header. if not create new entry.
		if aboveThreshold is False:
			listSeqs[seq1] = k.description
		if aboveThreshold is True:
			listSeqs[match] += "\t" + k.description

	noIdenticalOutfile = open("no_identicals." + output_file,"w")

	for k in listSeqs:
		noIdenticalOutfile.write(">" + listSeqs[k] + "\n" + k +"\n")

	print "\nIdentical sequences collapsed and written to %s" % (noIdenticalOutfile)

#function to return string difference
def similar(a, b):
	matching = sum(c==d for c, d in zip(a, b))
	totalLen = len(a)
	matching = float(matching)
	totalLen = float(totalLen)
	return matching/totalLen #return float of pairwise distance

#function to check that threshold is a value between 0 to 1.0
def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='''
Remove identical or similar sequences from a fasta alignment.

Input: Fasta

Requires: BLAST+, ACT and BioPython on your PATH
	''', formatter_class=RawTextHelpFormatter)

	#takes in the input files
	parser.add_argument('input', nargs="+", action="store", help="Specify at least 2 input files")

	parser.add_argument("-o", "--output", action="store", default=False, help="Output file. Default ...")
	parser.add_argument("-t", "--threshold",type=restricted_float, action="store", default="1.0", help="Maximum level  of sequence similarity to collapse 2 seqs'")


	#unused arguments
	parser.add_argument("-f", "--format", action="store", help="Specify output format e.g. fasta or phylip'")

	args = parser.parse_args()


	main() #run main script
