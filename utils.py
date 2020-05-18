#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

# to obtain a ref-read dataframe with count value for each sample
def get_ref_read_df(refset, count_table): #to obtain a ref-read dataframe with count value for each sample
    import pandas as pd
    id_dict = {ref.ID:{} for ref in refset}
    for refseq in refset:
        ref_id = refseq.ID
        read_ids = refseq.getReadIDs()
        id_dict[ref_id] = {read_id:count_table[read_id] \
                           for read_id in read_ids}
    return pd.DataFrame.from_dict(id_dict) #dictioary to a DF

# add reference`sequence strings to refseq objects
def fetch_refseq(RDPHOME, id_file_name, outfile_name, reffile_name):
    import subprocess
    import os
    subprocess.check_call(['java', '-jar', os.path.join(RDPHOME, \
                           'ReadSeq.jar'), 'select-seqs', id_file_name, \
                            outfile_name, 'fasta', 'Y', reffile_name],\
                            stdin=None, stdout=None, stderr=None, shell=False)
    return 1

def update_refseq(DF, reffile_name, refset):#add additional attributes to the refseq objects
    import string
    recs = open(reffile_name, 'r').read().strip('>').split('\n>')
    #Iterate all refseq objects for multiple tasks
    for rec in recs:
        lines = rec.split('\n')
        ID = lines[0].split()[0]
        seq = string.join(lines[1:], '').replace(' ', '')
        refset[ID].addSeq(seq)
        refset[ID].addReadCounts(DF)
        refset[ID].addRegs()
    return refset

## make a unique id for each consensus by adding a serial number to its template ID
class Name_proxy:
    counter = {}
    def get_assumed_id(self, ID):
        if not ID in self.__class__.counter:
            self.__class__.counter[ID] = 1
        else:
            self.__class__.counter[ID] += 1
        number = self.__class__.counter[ID]
        name = ID + '_' + str('0'*(3-len(str(number))) + str(number))
        return name

## classify consensus sequences in batch, fetch, and add assignments to Refseq objects
def classify_proxy(sample_id, RDPHOME, WD):
    import subprocess
    import os
    #classfy the consensus sequences
    subprocess.check_call(['java' , '-jar',\
            os.path.join(RDPHOME, 'classifier.jar'), \
            'classify',\
            '-f', 'fixrank', \
            '-o', os.path.join(WD, sample_id + '.cls'), \
            os.path.join(WD, sample_id + '_consensus.fasta')])
    return 1

def build_tree(seqfile_name, WD, RESDIR): #make a phylogenetic tree from template sequences
    import os
    aligned = os.path.join(WD, 'templates_mafft.fasta')
    tree = os.path.join(RESDIR, 'templates_mafft.tree')
    os.system('mafft --quiet --thread 4 FASTA > ALIGNED'.replace('FASTA', \
               seqfile_name).replace('ALIGNED', aligned))

    os.system('fasttree -quiet -nopr -nt -gtr < ALIGN > \
               TREE'.replace('ALIGN',aligned).replace('TREE', tree))
    return 1

def read_seq(seqfile_name):
    import string
    seq_dict = {}
    recs = open(seqfile_name, 'r').read().strip('>').split('\n>')
    for rec in recs:
        lines = rec.split('\n')
        ID = lines[0].split()[0]
        seq = string.join(lines[1:], '').replace(' ' , '')
        seq_dict[ID] = seq
    return seq_dict

def rev_complement(seq):
    anticodon = {"A":"T", "T":"A", "C":"G", "G":"C", "U":"A", \
                 "M":"K", "K":"M", "R":"Y", "Y":"R", "S":"S", \
                 "W":"W", "V":"B", "B":"V", "H":"D", "D":"H", \
                 "X":"X", "N":"N", "-":"-", ".":"."}
    revseq = ''
    while seq:#iterate over every character in the string
        base = seq[-1].upper()
        try:
            revseq = revseq + anticodon[base]
        except KeyError:
            revseq = revseq + base
        seq=seq[:-1]#delete the character positioned at[-1], which has been processed
    return revseq
