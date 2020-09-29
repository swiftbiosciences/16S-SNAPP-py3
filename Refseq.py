#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

#Class of templates to be used for consensus sequences
from operator import itemgetter
from itertools import groupby
import numpy as np
import string

class Refseq:
    def __init__(self, ID):
        self.ID = ID
        self.seq = '' #full length reference sequence
        self.reads = {}
        self.positions = []
        self.consensus = ''

    def addRead(self, read_info):#add read info: mapped positions
        ID = read_info['ID']
        start, end  = read_info['pos']
        positions = range(start, end) #mapped positions
        self.reads[ID] = positions #mapped positions
        for position in positions:
            if not position in self.positions:
                self.positions.append(position)
        self.positions.sort()

    def addPEread(self, ID, read_info_R1, read_info_R2):#add PE read info: two sets of positions
        s1, e1 = read_info_R1
        s2, e2 = read_info_R2
        positions = set(range(s1, e1)).union(set(range(s2, e2)))
        positions = list(positions)
        positions.sort()
        self.reads[ID] = positions #mapped positions
        for position in positions:
            if not position in self.positions:
                self.positions.append(position)
        self.positions.sort()

    def addPEread_bk(self, read_info_R1, read_info_R2):#add PE read info: two sets of positions
        ID = read_info_R1['ID'].split('_R')[0]
        s1, e1 = read_info_R1['pos']
        s2, e2 = read_info_R2['pos']
        positions = set(range(s1, e1)).union(set(range(s2, e2)))
        positions = list(positions)
        positions.sort()
        self.reads[ID] = positions #mapped positions
        for position in positions:
            if not position in self.positions:
                self.positions.append(position)
        self.positions.sort()

    def getReadIDs(self):
        return set(self.reads.keys())

    ## Execute the following operations when read mapping is complete and read count normalization is done
    def addSeq(self, seq): #the full length sequence of reference
        self.seq = seq.upper()

    def addRegs(self): #return continuously covered regions once the read mapping is complete
        self.baseRegs = [map(itemgetter(1), g) for k, g in groupby(enumerate(self.positions), lambda (i,x):i-x)]

    def addReadCounts(self, DF): #add minimized read counts attributable to this reference of this sample
        readIDs = self.reads.keys()
        self.count = DF.loc[readIDs, self.ID].to_dict() #count of reads by readID

        self.baseFreq = {pos:0 for pos in self.positions} #all mapped positions
        for readID, positions in self.reads.items():
            for pos in positions:
                self.baseFreq[pos] += self.count[readID]

    def addAssign(self, tax):#replace tax assignment with better resolved reads and eventually the emsemble's assignment
        self.tax = tax

    def getCountSum(self):#return the count sum of reads mapped to this reference
        return sum(self.count.values())

    def getMeanBaseCov(self):#return the average time a mapped base is covered by reads
        return np.array(self.baseFreq.values()).mean()

    def getAlignLen(self):
        length = 0
        for reg in self.baseRegs:
            length += len(reg)
        return length

    def getAlignPct(self):#fraction of all aligned positions 
        return round(len(self.positions)/float(len(self.seq)), 2)

    def getProxy(self):#return the proxy sequence formed by concatenating all mapped regions
        proxy = [self.seq[reg[0]:reg[-1]] for reg in self.baseRegs]
        return string.join(proxy, 'NNNNNNN')

if __name__ == '__main__':
    myref = Refseq(ID='r1', seq='GCAACTGGACTGGAA')
    myref.addReads(['asv1', 'asv2', 'asv3'])
    myref.addReads(['asv8'])
    myref.addBases(set([4, 5, 9, 3, 10]))
    myref.addBases(set([4, 5, 40, 59, 10]))
    print myref.ID
    print myref.seq
    print myref.readIDs
    print myref.getCount()
    print myref.bases
    print myref.getBaseCov()
    print myref.getRegs()
