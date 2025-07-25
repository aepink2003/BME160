#!/usr/bin/env python3
# Name: Ash Ediga (aediga)
# Group Members: None

import sequenceAnalysis
'''
    Simple unique subset finder class for tRna sequences

    getPowerSet() : find powerset of sequence
    getUnique() : find unique subsets by removing union of all powersets except this one
    getEssential() : compare all unique subsets and find minimum subsets
'''
class UniqueFinder:
    '''holds all power sets'''
    tRnaSets = []

    def __init__(self, sequence=""):
        self.sequence = sequence

    '''clean up input sequence'''
    def addSequence(self, newSequence):
        newSequence = newSequence.replace("_", "")
        self.sequence = newSequence.replace("-", "")

    '''get power set of the sequence'''
    def getPowerSet(self):
        powerset = []
        for i in range(len(self.sequence)):
            for j in range(i + 1, len(self.sequence) + 1):
                powerset.append(self.sequence[i:j])
        if powerset not in self.tRnaSets:
            self.tRnaSets.append(powerset)
        return powerset

    '''get all unique subsets of sequence'''
    def getUnique(self):
        allPowerSets = self.tRnaSets.copy()
        allPowerSets.remove(self.getPowerSet())
        uniqueSequence = set(self.getPowerSet()) - set().union(*allPowerSets)
        return uniqueSequence

    '''find the essential subsets in unique ones'''
    def getEssentials(self):
        nonessential = []
        uniqueList = list(self.getUnique())
        compareList = list(self.getUnique())

        for subset in uniqueList:
            compareList.remove(subset)
            for seq in compareList:
                if subset in seq:
                    if seq not in nonessential:
                        nonessential.append(seq)
            compareList.append(subset)

        return set(uniqueList) - set(nonessential)

'''find all indexes for essential subsets'''
def allindices(string, sub):
    i = 0
    listindex = []
    while i <= len(string):
        index = string.find(sub, i)
        if index == -1:
            break
        if index not in listindex:
            listindex.append(index)
        i += 1
    return listindex

'''save all powerset, find essential of each sequence, sort and print in correct format'''
def findUnique():
    fastaReader = sequenceAnalysis.FastAreader('bos-tRNA-7.fa')
    uniqueFinder = UniqueFinder()
    organizedFile = []

    for header, sequence in fastaReader.readFasta():
        organizedFile.append((header[8:11], header, sequence))
        uniqueFinder.addSequence(sequence)
        uniqueFinder.getPowerSet()
    organizedFile.sort()

    for aa, header, sequence in organizedFile:
        essentials = []
        uniqueFinder.addSequence(sequence)
        print(header)
        print(uniqueFinder.sequence)

        for essential in uniqueFinder.getEssentials():
            match = allindices(uniqueFinder.sequence, essential)
            for pos in match:
                tup = (essential, pos)
                essentials.append(tup)
        essentials.sort(key=lambda tup: tup[1])
        for sub in essentials:
            print("{0}{1}".format(sub[1]*".", sub[0]))

class main():
    findUnique()

main()

