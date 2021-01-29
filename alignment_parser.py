#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

#used to assemble K1 alignments from RDP pairwise alignment tool and return a
#sequence base dataframe assuming all base positions on the subject (1) as model
#positions and all inserts in the subject are non-model positions and are ignored

import pandas as pd
import sys
from itertools import islice

def remove_non_template_positions(template_base_list, read_base_list):
    read_base_list_mp = [] #base positions exist on the template
    for i in range(len(template_base_list)):
        if not template_base_list[i] == '-':
            read_base_list_mp.append(read_base_list[i])
    return read_base_list_mp

#read and convert each pair-wise alignments of all reads to the template
def get_align_array_bk(filename):
    align_list = []
    read_ids = []
    with open(filename, 'r') as infile:
        header = [line.strip() for line in islice(infile, 2)] #skip the first two lines
        while 1:
            try:
                lines = [line.strip() for line in islice(infile, 3)]
                read_id = lines[0].split('\t')[0]
                read_ids.append(read_id)
                read_align = list(lines[1][1:])
                template_align = list(lines[2][1:])
                read_align = remove_non_template_positions(template_align, read_align)

                align_list.append(read_align)
            except IndexError:
                break
    return align_list, read_ids

def get_align_array(alignment_lines):
    align_list = []
    read_ids = []
    i = 2
    line_num =  len(alignment_lines)
    while 1:
        read_id = alignment_lines[1].split('\t')[0]
        read_ids.append(read_id)
        read_align = list(alignment_lines[i + 1][1:])
        template_align = list(alignment_lines[i + 2][1:])
        read_align = remove_non_template_positions(template_align, read_align)
        align_list.append(read_align)
        i += 3
        if i + 3 > line_num:
            break
    return align_list, read_ids

#return seq_base_df
if __name__ == '__main__':
    array = get_align_array(sys.argv[1])
    print (array)
