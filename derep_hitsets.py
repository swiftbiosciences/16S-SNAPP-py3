#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

# the mini workflow to dereplidate blastn file prior to converge step
## simply dereplicate read sets by IDs to make a list of unique sets
## Ported to Python 3 on 20210106

import sys

#To make reference-to-read set hash keyed by reference seq IDs
def getSets(blastn):
    f = open(blastn, 'r')
    Hash = {}
    while 1:
        try:
            line = next(f)
            sID, qID = line.strip().split('\t')[0:2]
            if not sID in Hash:
                Hash[sID] = set([])
            Hash[sID].add(qID)
        except StopIteration:
            break
    return Hash

#Simply dereplicate the list of sets of read IDs
def getUniqSets(Hash):
    uniqSets = []
    dereped = {}
    for sID, qSet in Hash.items():
        if not qSet in uniqSets:
            uniqSets.append(qSet)
            dereped[sID] = qSet
    return dereped

def getDerepHitSets(blastn):
    sets = getSets(blastn)
    uniqSets = getUniqSets(sets)
    return uniqSets.keys()

if __name__ == '__main__':
    derepIDs = getDerepHitSets(sys.argv[1])
    for ID in derepIDs:
        print (ID)

