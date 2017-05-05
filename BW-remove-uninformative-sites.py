#!/usr/bin/env python

import sys, os
import argparse
from argparse import RawTextHelpFormatter
from collections import Counter #for counting characters in a string
from Bio import AlignIO

def main():    

    #print args.input[0] #Prints out input file name

    #define output file
    listRemovedColumns = []
    listFinalColumns = []
    
    #Check if output file given, if not use output prefix with input filename
    if not args.output:
        output_file = "output." + args.input[0]
    else:
        output_file = args.output
        
    alignment = AlignIO.read(args.input[0], "fasta")
    alignLen = alignment.get_alignment_length() #get length of alignment i.e. number of columns
    numTaxa = len(alignment) #count number of taxa
    
    print "\nProcessing an alignment with %s taxa and %s positions..." % (numTaxa, alignLen)
    
    
        #iterate across the length/columns of alignment
    for i in range(0,alignLen):
   
        #Print progress to stdout
        sys.stdout.write('\r')
        sys.stdout.write("Processing column %s" % (i))
        sys.stdout.flush()
   
        column = alignment[:, i] #reads in column number i as a string
        #print str(column) #print column in alignment
        
        dictColumnCounts = Counter(column) #Shows counts of caharacters
        dictColumnCountsValues = dictColumnCounts.values() 
        
        if (1 in dictColumnCountsValues) or (numTaxa in dictColumnCountsValues):
            listRemovedColumns.append(i) #List of rejected columns
            
            
        else:
            listFinalColumns.append(i) #List of accepted columns
       
    
    sys.stdout.flush()

    
    finalAlignment = alignment[:, listFinalColumns[0]:listFinalColumns[0]+1] #initialise new Multiple seq alignment with the first column from listFinalColumns
    #print listFinalColumns[0]
    
    for j in listFinalColumns[1:]: #iterate through list of accepted columns, skipping the first one which is already used to initialise the object
        
        #print j
        
        #Append the columns to the final alignment
        finalAlignment = finalAlignment + alignment[:, j:j+1]
        
        #Print progress to stdout
        sys.stdout.flush()
        sys.stdout.write('\r')
        sys.stdout.write("Appending column %s" % (j))
        sys.stdout.flush()
    
    #print finalAlignment #Print alignment object
    print "\nFinal alignment contains %s columns" % (finalAlignment.get_alignment_length())
            
        
    AlignIO.write(finalAlignment, output_file, "fasta")
    
    print "\nRemoved a total of %s columns, with %s columns informative columns remaining" % (len(listRemovedColumns), len(listFinalColumns))
        
                
#subroutine to count occurrences of characters in a string
def chunk_string(s, n):
    return [s[i:i+n] for i in range(len(s)-n+1)]




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
Remove uninformative sites from fasta alignment.
    
Input: Fasta
    
Requires: BLAST+, ACT and BioPython on your PATH
    ''', formatter_class=RawTextHelpFormatter)
        
    
    #takes in the input files
    parser.add_argument('input', nargs="+", action="store", help="Specify at least 2 input files")
    
    parser.add_argument("-o", "--output", action="store", default=False, help="Output file. Default ...")

    #unused arguments
    parser.add_argument("-f", "--format", action="store", help="Specify output format e.g. fasta or phylip'")



    main() #run main script

