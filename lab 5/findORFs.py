#!/usr/bin/env python3
# Name: Ash Ediga (CATS account username)
# Group Members: “None”
#design 
#ORFfinder class is a part of the sequenceAnalysis module 
#init
#reads in the DNA sequence and creates a list of ORFs found in the sequence
#findORF
#searches for codons by dividing sequence into 3 nucleotide slices using frames 
#if start codon is found the location is stored in a list 
#if stop is found before start, position is changed to 1 
#if both codons are found indicating start & stop the ending position sets to +3 on a frameshift
# if start codon with stop if found the position is set to end of the input sequence
#reverse start has the reverse done to it 
#reverseComplement: creates the complement using an opposite sequence
#main
#reads using commandline and the sequence from fasta file 
#prints in requested format
from sequenceAnalysis import OrfFinder, FastAreader
########################################################################
# CommandLine
########################################################################
class CommandLine():
    
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (100,200,300,500,1000), default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', 
                                 help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Main
# Here is the main program
# 
#
########################################################################
def main(inCL=None):
    """Reads in a fasta file and outputs the ORFs frame, start, stop, and length position on a output file."""
    if inCL is None:
        thisCommandLine = CommandLine()
        if thisCommandLine.args.longestGene:
            fastaFile = FastAreader()
            for header, sequence in fastaFile.readFasta():
                print(header)
                orfData = OrfFinder(sequence)
                orfData.findReverseOrfs()
                orfData.findOrfs()
                filteredList = filter(lambda orf: orf[3] > myCommandLine.args.minGene, orfData.orfs)  # Filters out the ORFS based on the minGene arg.
                for frame, start, stop, length in sorted(filteredList, key=lambda orf: orf[3], reverse = True):  # Sort list of ORFs by length.
                    print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(frame, start, stop, length))
    else:
        thisCommandLine = CommandLine(inCL)
    print(thisCommandLine.args)

if __name__ == "__main__":
    main()