#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

import os
import string
import concurrent.futures
import alignment_parser

def align_seqs(refseq, pe_seq_dict, WD):
    read_ids = refseq.reads.keys()
    read_ids.sort()
    read_seqs = [pe_seq_dict[read_id] for read_id in read_ids]
    read_count_list = []
    #read count for both R1 and R2 reads
    for read_id in read_ids:
        read_count_list.append(refseq.count[read_id])
        read_count_list.append(refseq.count[read_id])
    template_id = refseq.ID
    template_seq = refseq.seq
    read_filename = os.path.join(WD , template_id + '_reads.fasta')
    template_filename = os.path.join(WD, template_id + '_temp.fasta')

    #write the read sequences to a fasta file
    with open(read_filename, 'w') as out:
        for i in range(len(read_ids)):
            r1_seq, r2_seq = read_seqs[i].split('N'*10)

            #The PE reads have to be aligned separately
            out.write('>' + read_ids[i] + '_R1' + '\n' + r1_seq + '\n')
            out.write('>' + read_ids[i] + '_R2' + '\n' + r2_seq + '\n')
    #write the template sequence to a fasta file
    with open(template_filename, 'w') as out:
        out.write('>' + template_id + '\n' + template_seq)

    aligned = os.path.join(WD, template_id + '_rdpaligned.txt')

    #Run alignment tool align each read to the template in glocal mode
    #collect alignment from stdout
    aligned = os.popen('java -jar ~/RDPTools/AlignmentTools.jar pairwise-knn \
               -k 1 READS TEMPLATE'.replace('READS', read_filename).\
               replace('TEMPLATE', template_filename)).readlines()

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
    consensus = string.join(consensus, '')

    #Replace common gaps with a fixed length (7) of 'N's
    consensus = string.join([region for region in consensus.split('-') \
            if not region == ''], 'N'*7)
    return consensus

def consensus_loader(refseq, pe_seq_dict, WD):
    alignment, count_series = align_seqs(refseq, pe_seq_dict, WD)
    consensus = get_consensus(alignment, count_series)
    seq = '>' + refseq.ID + '\n' + consensus.strip().strip('N')
    #refseq.consensus = consensus #load consensus sequences
    template_id = refseq.ID
    tmp_file_prefix = os.path.join(WD, template_id)
    os.system('rm %s*'%tmp_file_prefix) #remove the temperoray files
    return seq

def load_consensus(refset, pe_seq_dict, WD):
    all_consensus = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        results = [executor.submit(consensus_loader, refseq, pe_seq_dict, WD) \
                   for refseq in refset.values()]
        for f in concurrent.futures.as_completed(results):
            all_consensus.append(f.result())
    return string.join(all_consensus, '\n')
