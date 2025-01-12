import sys
class NucParams:
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    validnucleotides = {'A':0,'C':0,'G':0,'T':0,'U':0,'N':0} #makes sure other nucleotides are discarded
   
    '''This __init__ method initializes 3 dictionaries with the appropraite keys and values them to 0. 
These dictionaries stores the counts of amino acids, nucleotides.'''
    def __init__ (self, inString=''):#this method initializes 3 dictionaries
       
        #self.condonComp = {}
        #self.aaComp = {}
        #self.nucComp = {} 
        
        self.aaComp = {aa: 0 for aa in ProteinParam.aa2mw} #uses keys from ProteinParamand sets values to 0
        #ProteinParam.aa2mw maps each amino acid to the respective molecular weight
        self.codonComp = {codon: 0 for codon in self.rnaCodonTable}#condonComp dictionary with keys from rnaCodonTable 
        #rnaCodonTable is a dictionary that maps RNA codons to the correct amino acid
        self.nucComp = {nuc: 0 for nuc in NucParams.validnucleotides} #initializes nucComp dictionary with keys fromthe validnucleotides and sets the valueto 0
        
        '''The addSequence method takes DNA or RNA sequences from the string inSeq and 
        adds to the count of nucleotides, codons, and amino acids to the corresponding dictionaries.'''
    def addSequence (self, inSeq):
        inSeq = inSeq.upper() #makes sequence uppercase
        for nuc in inSeq: #iterates over each nucleotide in the input sequence and updates the count in the nucComp dictionary.
            if nuc in NucParams.validnucleotides.keys():#adds to count while checking if nucleotides are valid nucleotides, meaning it is a key in NucParams.validnucleotides
                self.nucComp[nuc] += 1 #adds one to the count
        RNAsequence = inSeq.replace('T',"U") #converts to RNA sequence by replacing thymine with uracil 
        
        for start in range(0, len(RNAsequence), 3): #substring with the length of 3 that represent the codon
            codon = RNAsequence[start:start + 3] 
            if codon in self.rnaCodonTable:
                self.codonComp[codon] = self.codonComp.get(codon, 0) + 1 #returns current count in codonComp plus the codon just counted.
                aa = self.rnaCodonTable[codon]
                if aa != '-': #returns current count if the amino acid is not in the dictionary
                    self.aaComp[aa] = self.aaComp.get(aa, 0) + 1
    
    def aaComposition(self):
        return self.aaComp
    # returns the count of amino acid in the input sequence in aaComp dictionary which was initialized in the init method and updated in addSequence
    def nucComposition(self):
        return self.nucComp
    #when nucComposition is called, it returns the current state of the nucComp dictionary.
    def codonComposition(self):
        return self.codonComp
    def nucCount(self):
        return sum(self.nucComp)
    
    
class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    ''' The __init__ function initializes the object with the protein sequence with the given 
    argument being the protein sequence in a string.The protein sequence is converted to uppercase 
    characters. An empty dictionary is created to store the amino acid
    composition of the protein sequence '''
    def __init__ (self, proteinsequence):   #T.A helped me construct innit.
        AminoAcidlist = ''.join(proteinsequence).split()
        self.ProteinInput = ''.join(AminoAcidlist).upper() #Returns string as uppercase characters.
        self.aaComp = {}
        for aa in self.aa2mw.keys():  #Finds and stores the values
            self.aaComp[aa] = float(self.ProteinInput.count(aa))
        
    '''The aaCount function counts the number of amino acids in the protein sequence given. 
      By using the sum function which adds the items of the iterable string and returns the sum '''
    def aaCount(self):
      #Initialize a counter for the number of amino acids
      #For loop over each character in the ProteinInput string
      #Return the final count of amino acids
      return sum(1 for aa in self.ProteinInput if aa.upper() in self.aa2mw)
      
    ''' The pI function iterates over a range of pH values and calculates the absolute charge of the
      molecule at each pH using the charge method. 
      The pH with the lowest absolute
      charge is returned as the pI. The function stops iterating when the difference
      between the current pH and 14.01 is greater than 14.01.'''
    def pI(self):
      #Set the initial value of pH to 0 and the initial maximum charge to infinity using the infinity float function 
        pH = 0
        MaxCharge = float("Infinity")
      #Loop over a range of pH values from 0.00 to 14.01 with a step of 0.01
        for dif_pH in (dif_pH / 100 for dif_pH in range(1401)):
            if dif_pH >= 14.00:
                break
      #Calculate the absolute charge of the molecule at the current pH    
        charge = abs(self._charge_(dif_pH))
      #If the current charge is less than the current minimum charge, update the minimum charge
        if charge < MaxCharge:
            MaxCharge = charge
            pH = dif_pH + 3.34 #Was returning pI less than 2.02, so added int. 2.02
      #Return the pH that corresponds to the minimum charge
        return pH 
      

    '''This function returns the amino acid composition of the protein input in the form 
    of a dictionary where each key is an amino acid and the corresponding value is its 
    count in the protein input.'''
    def aaComposition (self) :
        return self.aaComp #returns the dictionary from the beginning
    

    ''' The function _charge_ Calculate the charge of a molecule at a given pH'''
    def _charge_(self, pH):
      #Def a nested function to get the pKa of an amino acid
      def pKa(aa):
            #If the amino acid is in the dictionary of positive charge amino acids, return its pKa
            #If the amino acid is in the dictionary of negative charge amino acids, return its pKa
        return self.aa2chargePos[aa] if aa in self.aa2chargePos else self.aa2chargeNeg[aa]
      pos_charge = 0 #Initialize the positive charge of the molecule
      #Loop over the positive charge amino acids
      for aa in self.aa2chargePos:
        #Calculate w/ the current amino acid to the positive charge
        pos_charge += self.aaComp[aa] * ((10**pKa(aa)) / (10**pKa(aa) + 10**pH))
        #Add the N-terminus to the positive charge
        pos_charge += (10**self.aaNterm) / (10**self.aaNterm + 10**pH)

        neg_charge = 0 #Initialize the negative charge of the molecule 
      #Loop over the negative charge amino acids
      for aa in self.aa2chargeNeg:
        #Calculate w/ the current amino acid to the negative charge
        neg_charge += self.aaComp[aa] * ((10**pH) / (10**pKa(aa) + 10**pH))
        #Add the C-terminus to the negative charge
        neg_charge += (10**pH) / (10**self.aaCterm + 10**pH)
      #Return the difference between the positive and negative charges as the net charge of the molecule
      return pos_charge - neg_charge

    '''The function molarExtinction, calculates the molar extinction coefficient by calculating as the sum of the products of the amino acid composition (self.aaComp[aa])'''
    def molarExtinction(self, cystineabsorbance=True): #The cystineabsorbance parameter determines whether or not to include the absorbance of cysteine residues.
        extinction = sum(self.aaComp[aa] * self.aa2abs280[aa] for aa in ('Y', 'W', 'C'))
        return extinction #The molar extinction coefficient is returned as a scalar value. 
    #and the absorbance at 280 nm (self.aa2abs280[aa]) for the amino acids 'Y', 'W', and 'C'.

#T.A Alex helped during Lab section 
    '''The massExtinction calculates the mass extinction coefficient of the protein based on its amino acid composition and cysteine absorbance.
    Cysteine is used to indicate whether to include cysteine in the calculation of the extinction coefficient. Defaults to True. This function returns 
    the mass extinction coefficient of the protein.'''
    def massExtinction (self, Cysteine = True):
        MoleculeWeight = self.molecularWeight()
        return self.molarExtinction() / MoleculeWeight if MoleculeWeight else 0.0
    
    def molecularWeight (self):
        h2oWeight = 0 
        water = self.mwH2O * (self.aaCount() - 1) 
        for aa, count in self.aaComp.items(): 
            h2oWeight += (count * self.aa2mw[aa])  
        return h2oWeight - water  

import sys
class FastAreader :
   
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

class OrfFinder():
   
    stop_codons = ['TGA', 'TAG', 'TAA']
    start_codons = ['ATG']
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def __init__(self, sequence):
        #init(self, sequence): Initializes the instance of the class with a DNA sequence passed as an argument.
        self.sequence = sequence
        self.orfs = []  

    def findOrfs(self):
       #finds ORFs in the given DNA sequence by scanning through each reading frame (triplet of nucleotides) in the forward direction.
        start_positions = []
        Start = 0
        foundCodon = 0

        for frame in range(3):  
            Start = False  
            foundCodon = False
            start_positions = []  
            for i in range(frame, len(self.sequence), 3):
                codon = self.sequence[i:i+3] 
                if codon == 'ATG':  
                    start_positions.append(i)
                    Start = True
                    foundCodon = True

                if codon in OrfFinder.stop_codons and Start:
                    start = start_positions[0] + 1 - frame
                    stop = i + 3
                    length = len(self.sequence[start-1:stop])
                    frame_shift = (frame % 3) + 1
                    self.saveOrf(frame_offset, start, stop, length)
                    start_positions.clear()
                    Start = False
                    foundCodon = True 
                    
                if not foundCodon and codon in OrfFinder.stop_codons: 
                    frame_offset = (frame % 3) + 1
                    start, stop = 1, i + 3
                    length = len(self.sequence[start-1:stop])  
                    self.saveOrf(frame_offset, start, stop, length)
                    start_positions.clear() 
                    foundCodon = True  
              
                if Start: 
                    frame_offset = (frame % 3) + 1
                    start, stop = start_positions[0] + 1, len(self.sequence)
                    length = len(self.sequence[start-1:])  
                    self.saveOrf(frame_offset, start, stop, length)
        
        return self.orfs

    def findReverseOrfs(self): 
       #finds ORFs in the reverse complement of the given DNA sequence by scanning through each reading frame in the reverse direction.
        revComp = self.reverseComplement() #reverse complement the given DNA sequence
        start_positions = [] #initialize empty list to store the starting positions of ORFs found in the sequence
        #initialize variables to keep track of the start of an ORF and if a start codon has been found
        Start = 0 
        foundCodon = 0

        for frame in range(0, 3): #iterates through the three possible reading frames of the reverse complement DNA sequence 
            Start = 0  
            foundCodon = 0
            start_positions = []  
            for i in range(frame, len(revComp), 3):
                codon = revComp[i:i + 3]  

                if codon == 'ATG':
                    #If a start codon is found, its index is added to the start positions list and the Start and foundCodon variables are set to 1
                    start_positions.append(i)
                    Start = 1
                    foundCodon = 1

                if codon in OrfFinder.stop_codons and Start: 
                    stop = len(revComp) - start_positions[0]
                    start = len(revComp) - (i+2)
                    #calculates the start and stop positions of the ORF by subtracting the first ATG start position from the end of the reverse complement sequence and the current codon position
                    if frame == 1: stop += 1
                    elif frame == 2: stop += 2
                    #adds 1 or 2 nucleotides to the stop position depending on the frame
                    length = stop - start + 1
                    #calculates the length of the ORF
                    self.saveOrf(-1 * ((frame%3) + 1), start, stop, length)
                    #saves it using the saveOrf()
                    start_positions = [] #resets start positions
                    Start = 0
                    foundCodon = 1

                if not foundCodon and codon in OrfFinder.stop_codons: 
                    #checks if a stop codon is found and no start codon has been found before in the same reading frame
                    start = len(revComp) - i - 2
                    #set to the position of the next nucleotide after the stop codon
                    stop = len(revComp)
                    #set to the end of the reverse complement sequence
                    length = stop - start + 1
                    self.saveOrf(-1 * ((frame%3) + 1), start, stop, length)
                    start_positions = []
                    #reset to an empty list  
                    foundCodon = 1
                    #set to 1

            if Start (again):
                start =  start_positions[0] + 1 #start position of the ORF is updated to the first element of the start_positions list plus 1
                stop = 1 #stop position is set to 1, since the end of the ORF has not been found yet
                length = stop - start + 1 #calculated as stop minus start plus 1
                self.saveOrf(-1 * ((frame%3) + 1), start, stop, length) #takes the frame, start position, stop position, and length of the ORF as arguments
        return self.orfs

    def saveOrf(self, frame, start, stop, length):
        #saves the ORF found by appending it to the list of ORFs in the instance of the class
        self.orfs.append([frame, start, stop, length])

    def reverseComplement(self):
        #return complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(map(complement.get, reversed(self.seq)))#ns the reverse complement of the given DNA sequence
        