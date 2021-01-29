#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

import os
import time
import concurrent.futures
import alignment_parser
import numpy as np
RDPHOME = os.environ['RDPHOME']

def align_seqs(refseq, pe_seq_dict, tmp_dir):
    read_ids = refseq.reads.keys()
    #read_ids.sort() # changed as the next line for Python 3
    read_ids = sorted(read_ids)
    read_seqs = [pe_seq_dict[read_id] for read_id in read_ids]
    read_count_list = []
    #read count for both R1 and R2 reads
    for read_id in read_ids:
        read_count_list.append(refseq.count[read_id])
        read_count_list.append(refseq.count[read_id])
    template_id = refseq.ID
    template_seq = refseq.seq
    read_filename = os.path.join(tmp_dir , template_id + '_reads.fasta')
    template_filename = os.path.join(tmp_dir, template_id + '_temp.fasta')

    #write the read sequences to a fasta file
    out1 = open(read_filename, 'w')
    for i in range(len(read_ids)):
        r1_seq, r2_seq = read_seqs[i].split('N'*10)

        #The PE reads have to be aligned separately
        out1.write('>' + read_ids[i] + '_R1' + '\n' + r1_seq + '\n')
        out1.write('>' + read_ids[i] + '_R2' + '\n' + r2_seq + '\n')
    out1.close()

    #write the template sequence to a fasta file
    out2 = open(template_filename, 'w')
    out2.write('>' + template_id + '\n' + template_seq)
    out2.close()

    #Run alignment tool align each read to the template in glocal mode
    #collect alignment from stdout
    aligned = os.popen('java -jar RDPHOME/AlignmentTools.jar pairwise-knn \
               -k 1 READS TEMPLATE'.replace('READS', read_filename).\
               replace('TEMPLATE', template_filename).replace('RDPHOME', RDPHOME)).readlines()
    return (aligned, read_count_list)

def get_consensus(rdp_k1_alignment, count_list):
    import pandas as pd
    #convert the alignment into a dataframe of sequence/bases and frequencies 
    #of all reads into a frequence matrix exactly matching the dataframe of 
    #sequence/bases

    count_array = [] #list of counts per sequence (row), and base (column)
    aligned_array, read_ids = alignment_parser.get_align_array(rdp_k1_alignment)
    for seq, count in zip(aligned_array, count_list):
        base_count = []
        for base in seq:
            if base in ['A', 'T', 'G', 'C']:
                base_count.append(count)
            else:
                base_count.append(0)
        count_array.append(base_count)
    aligned_df = pd.DataFrame(aligned_array)
    count_df = pd.DataFrame(count_array)

    #list of consensus bases
    consensus = []
    for i in range(len(aligned_df.columns)):
        pos_df = pd.concat([aligned_df.iloc[:, i], count_df.iloc[:, i]], \
                axis=1, keys=['base', 'count'])

        ##the base that has the highest frequency
        base = pos_df.groupby('base').sum().idxmax()['count']
        consensus.append(str(base))
    consensus = ''.join(consensus)

    #Replace common gaps with a fixed length (7) of 'N's
    consensus = "NNNNNNN".join([region for region in consensus.split('-') \
            if not region == ''])
    return consensus

def consensus_loader(sample_id, refseq, pe_seq_dict, tmp_dir):
    alignment, count_series = align_seqs(refseq, pe_seq_dict, tmp_dir)
    size = round(np.array(count_series).sum(), 3) # total read count of this consensus sequence
    consensus = get_consensus(alignment, count_series)
    seq = '>' + refseq.ID + ';' + 'sample_id=' + sample_id + ';size=%s'%size + '\n' + consensus.strip().strip('N')
    return seq

def load_consensus(sample_id, refset, pe_seq_dict, wd):
    all_consensus = []
    tmp_dir = os.path.join(wd, sample_id.split('_R1')[0])
    os.system('mkdir %s'%tmp_dir)
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        results = [executor.submit(consensus_loader, sample_id, refseq, pe_seq_dict, tmp_dir) \
                   for refseq in refset.values()]
        for f in concurrent.futures.as_completed(results):
            all_consensus.append(f.result())
    os.system('rm -fr %s'%tmp_dir) #remove the tmp directory
    return "\n".join(all_consensus)
