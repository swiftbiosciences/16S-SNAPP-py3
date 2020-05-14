#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

import sys
import string

#extract ASV sequences to a fasta file and replace the sequence strings with IDs
def split_asv_count(tbl, prefix):
    f = open(tbl, 'r')
    tblout = open(prefix + '_count.csv', 'w')
    single = open(prefix + '_seq.fasta', 'w')
    pair = open(prefix + '_PE.fasta', 'w')
    count = 0
    tblout.write(f.next().replace('"', ''))
    while 1:
        try:
            line = f.next()
            count += 1
            cols = line.strip().split(',')
            PE = cols[0].replace('"', '')
            R1R2 = PE.split('NNNNNNNNNN')
            ID = 'asv_%s'%count
            if len(R1R2) == 2:
                R1, R2 = R1R2

                fasta1 = '>' + ID + '_R1' + '\n' + R1 + '\n'
                fasta2 = '>' + ID + '_R2' + '\n' + R2 + '\n'
                single.write(fasta1 + fasta2)
                fasta = '>' + ID + '\n' + PE
                pair.write(fasta + '\n')
                cols[0] = ID
                tblout.write(string.join(cols, ',') + '\n')
            else:
                print 'error'
                sys.exit()

        except StopIteration:
            break
    single.close()
    pair.close()
    tblout.close()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'get_asv_files.py asvTable outprefix'
        sys.exit()
    asvTable = sys.argv[1]
    outprefix = sys.argv[2]
    split_asv_count(asvTable, outprefix)
