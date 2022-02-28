#!/usr/bin/env python
#Used to split fasta file into multiple small files of specified size

import sys
import os
if len(sys.argv) != 4:
    print 'Usage: $splitFasta.py fastaFile numberSeqInEach outdir'
    sys.exit()

f = open(sys.argv[1], 'r')
size = int(sys.argv[2])
outdir = os.path.abspath(sys.argv[3]) #output directory

if not os.path.isdir(outdir):
    os.system('mkdir %s'%outdir)

count = 0
File = 1
stem = os.path.join(outdir, os.path.basename(sys.argv[1]).split('.')[0])
out = open(stem + '_out%s.fas'%File, 'w')#output stream
line = f.readline()
while 1:
    if len(line.strip()) == 0:
        line = f.readline()
    if line.strip()[0] == '>':
        count += 1
        if count > size:
            out.close()
            File += 1 #next file
            out = open(stem + '_out%s.fas'%File, 'w')#output stream
            count = 1
    out.write(line)
    line = f.readline()
    if not line:
        break
out.close()


