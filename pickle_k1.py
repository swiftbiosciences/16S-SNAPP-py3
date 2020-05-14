#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

#prepare rdp seqmatch results hits as the dictionary and save as pickle for later use
import sys
import pickle

def pickle_k1(K1):
    Dict = {}
    for line in open(K1, 'r').readlines():
        asvID, templateID = line.strip().split('\t')[:2]
        Dict[asvID] = templateID
        f = open('asv_PE_K1.pkl', "wb")
        pickle.dump(Dict, f)
if __name__ == '__main__':
    if not len(sys.argv) == 2:
        print 'pickleK1.py asv_PE_K1.txt'
        sys.exit()
    pickle_k1(sys.argv[1])
