#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

#Take the blast results and coverge reads to their likely templates
## Ported to Python 3 on 20210108

import sys
import os
import pandas as pd
from blastn_parser import *
from utils import *
from minimize_var import *
from load_consensus import *
from lineage_parser import *
import concurrent.futures
import timeit

#Process each sample separately: each sample is represented by a read (row)-
#refence (column) count dataframe
def assign_read_counts(sampleID, readCounts, hits_info):
    print ('Allocating multi-mapped read counts for sample', sampleID)
    print ('Total mapped references:', len(hits_info))
    refSet_start_time = timeit.default_timer()

    ##1. instantiate refseq objects to associate reads
    refset, refseq_PE_list = initiate_refset(hits_info)
    print ('Pre-reduction refset size:', len(refset))
    print ('Number of matched PEs:', len(refseq_PE_list))
    print (sampleID, 'refset total time', timeit.default_timer() - refSet_start_time)

    ##2. Converge the refset to a much smaller set containing the greatest gene 
    ##   coverage and account for reference-mapped read counts
    refset = converge_ref(refset)
    print ('Minimized refset size:', len(refset))

    ##3. make a DataFrame of ref-read and count
    DF = get_ref_read_df(refset, readCounts).fillna(0)

    ##4. Sort row (reads) of DF by the nonzero read (max in this case) count of each rows
    #DF = DF.ix[DF.max(axis=1).sort_values(ascending=True).index,:] #changed to next line for Python 3
    DF = DF.loc[DF.max(axis=1).sort_values(ascending=True).index,:]

    ##5. Make a deep copy of this DF to hold assigned read counts for this sample
    DF_assigned = DF.copy(deep=True)
    print ('The sample ASV-reference dataframe', DF.shape)

    ##6. slice a subset (rows) where each read is mapped to more than one reference
    readIDs_tbd = DF.loc[(DF.astype(bool).sum(axis=1) > 1) , :].index

    ##7. Set all multi-mapped ref-read count to small number 0.001 They will be filled later with
    #the optimized values through minimization
    DF_assigned.loc[(DF_assigned.astype(bool).sum(axis=1) > 1) , :] = 0.001
    print ('There are:', len(readIDs_tbd), 'multi-mapped asv representatives')

    ##8. Iterate each multi-mapped read to optimize its count allocation to each
    ##   mapped reference by minimizing the sum of variance across all references 
    ##   it is mapped to. It is the most time consuming step and has to be done
    ##   sequencially
    for readID in readIDs_tbd: #iterate over each multi-mapped read
        df = DF.loc[:, (DF.loc[readID] > 0)] #dataframe slice containing all reference columns that have nonzero values with this read
        df = df.loc[(df != 0).any(1), :] #drop read rows that contain all zeros
        df_assigned = DF_assigned.loc[df.index, df.columns] #make a same slice the dataframe to hold normalized values

        rIDs = df.index.to_list() #all readIDs in this matrix (index)
        mask = [-1 for ID in rIDs]
        mask[rIDs.index(readID)] = df.loc[readID, :][0]

        # submit the dataframe and mask for allocating read counts by minimizing
        # the sum of variance of read counts across all reference mapped to this read
        df_assigned = minimize_var(df_assigned.T, mask)

        # update read count to the assigned. Do we need this step? Should assignments to the setset of DF be inplace?
        DF_assigned.update(df_assigned.T)

    ##9. Change post-minimization refset from a list to a dictionary keyed by their IDs
    refset = {refset[i].ID:refset[i] for i in range(len(refset))}

    ##10. Fetch refseq strings
    idList = os.path.join(WD, sampleID + '_idlist.txt')
    out = open(idList, 'w')
    out.write("\n".join(refset.keys()))
    out.close()

    seqFile = os.path.join(WD, sampleID + '_ref.fasta')
    fetch_refseq(RDPHOME, idList, seqFile, ref_db_path)
    refset = update_refseq(DF_assigned, seqFile, refset)

    return refset, refseq_PE_list

## Obtain alignment length distributions: per refseq and normlized by readCounts
def get_len_stats(sample_id, refset):
    length_per_refseq = []
    length_norm_readcount = []
    for refID in refset:
        length = refset[refID].getAlignLen()
        readCounts = int(refset[refID].getCountSum())
        length_per_refseq.append(length)
        length_norm_readcount.extend([length]*readCounts)
    length1 = '\n    ' + pd.Series(length_per_refseq).describe().to_string().replace('\n', '\n    ')
    length2 = '\n    ' + pd.Series(length_norm_readcount).describe().to_string().replace('\n', '\n    ')
    with open(log, 'a') as logfile:
        logfile.write('\n    Alignment length stats per template in ' + sample_id )
        logfile.write(length1 + '\n')
        logfile.write('\n    Alignment length stats per asv in ' + sample_id )
        logfile.write(length2 + '\n')
    return pd.Series(length_per_refseq).describe(), pd.Series(length_norm_readcount).describe()

##cluster consensus and asv sequences from all samples and make an abundance table with taxonomic assignments
def get_cluster_tbl(wd):
    all_seqs = os.path.join(wd, 'consensus_asv.fasta')
    os.system('cat %s/*_consensus.fasta > %s'%(wd, all_seqs))
    uc = os.path.join(wd, 'all.uc')
    reps = os.path.join(wd, 'reps.fasta')
    cmd = 'vsearch --cluster_size %s --strand both \
            --id 0.97 --uc %s --centroid %s'%(all_seqs, uc, reps)
    os.system(cmd)

def combine_lineage_count(sample_id, refset, unmapped_list):
    refseqs = refset.values()
    Hash = {} #to hold the feature counts collapsed by lineage

    feature_counts_dict = {} # the count series of associated and unassociated reads using the IDs of the templates or K1 nearest reference
    feature_taxa_dict = {} # the taxonomic assignments of the features
    feature_template_seqs_dict = {} # the template sequence strings

    #first fetch lineage-count from refseq objects
    #nameproxy = Name_proxy() #instantiate consensus naming object
    for refseq in refseqs:
        ID = refseq.ID
        lineage = refseq.tax.strip()
        count = refseq.getCountSum()
        template_seq = refseq.seq
        template_id = ID + ';sample_id=' + sample_id  #use modified template ID to avoid redundancy with other samples
        #template_id = nameproxy.get_assumed_id(ID) #use modified template ID to avoid redundancy with other samples
        feature_counts_dict[template_id] = count
        feature_taxa_dict[template_id] = lineage
        feature_template_seqs_dict[template_id] = template_seq
        if not lineage in Hash:
            Hash[lineage] = count
        else:
            Hash[lineage] += count

    #second fetch lineage-counts of unmapped asvs
    for asvID in unmapped_list:
        lineage = pe_lineage_dict[asvID].strip()
        count = rDF.loc[asvID, sample_id]
        try:
            k1_ID = K1_dict[asvID]
        except KeyError:
        #this read will not be included in the final tables because it is
        #beyond the seqmatch cutoff

            continue
        #Assign count and lineage of asvs to reference IDs
        if not k1_ID in feature_counts_dict.keys():
            feature_counts_dict[k1_ID] = count
            feature_taxa_dict[k1_ID] = lineage
        else:
            if lineage.count(';') > feature_taxa_dict[k1_ID].count(','): #pick the better assignmentSample
                feature_counts_dict[k1_ID] += count
                feature_taxa_dict[k1_ID] = lineage #replace with this better asssignment
            else:
                feature_counts_dict[k1_ID] += count
        if not lineage in Hash:
            Hash[lineage] = count
        else:
            Hash[lineage] += count
    return (Hash, feature_counts_dict, feature_taxa_dict, feature_template_seqs_dict)

def converge(sample_id): #converge each sample calculating all required attributes of the refset
    sample_count = rDF.loc[:, sample_id] #asv count series of this sample

    sample_count_series = sample_count[sample_count > 0]#drop 0's count series of this sample
    sample_read_count = len(sample_count_series)

    #load previously pickled blastn results as a dictionary for this sample
    blastn_pickle = os.path.join(WD, 'pickle/' + sample_id + '.pkl')
    blastn_match_dict = pd.read_pickle(blastn_pickle)

    #call assign_read_counts function to converge the asv PEs
    refset, mapped_ASV_IDs = assign_read_counts(sample_id, sample_count_series, blastn_match_dict)

    unmapped_ASV_IDs = list(sample_count_series.index.drop(mapped_ASV_IDs))
    mapped_ASV_count =len(mapped_ASV_IDs)
    with open(log, 'a') as logfile:
        logfile.write('    Count of total unique ASVs in ' + sample_id + ': ' + \
                       str(sample_read_count) + '\n')
        logfile.write('    Count of mapped unique ASVs in ' + sample_id + ': ' + \
                       str(mapped_ASV_count) + '\n')
        logfile.write('    Percentage of mapped unique ASVs in ' + sample_id + \
                      ': ' + str(round(mapped_ASV_count/float(sample_read_count)*100, 2)) + '\n')

    consensus_seqs = load_consensus(sample_id, refset, pe_seq_dict, WD)
    seq_name = os.path.join(WD, sample_id + '_consensus.fasta')
    with open(seq_name, 'w') as conseq: #write consensus seqs to a file
        conseq.write(consensus_seqs)
    classify_proxy(sample_id, RDPHOME, WD) #classify consensus sequences
    refset = add_lineages(WD, sample_id, pe_lineage_dict, refset)#load the best assignment from asv_PEs and the proxies

    #Append unmapped ASVs to the consensus sequence file with size information
    with open(seq_name, 'a') as asv:
        asv.write('\n')
        for asv_id in unmapped_ASV_IDs:
            size =  rDF.loc[asv_id , sample_id]
            asv.write('>' + asv_id + ';sample_id=' + sample_id + ';size=%s'%size + '\n' + pe_seq_dict[asv_id] + '\n')

    get_len_stats(sample_id, refset) #alignment length distributions of mapped reads
    collapsed_series, \
    feat_count_series, \
    feat_tax_series, \
    feat_tempSeq_dict \
    = combine_lineage_count(sample_id, refset, list(unmapped_ASV_IDs))#get abundance dictionary of all lineages in this sample

    return  (pd.DataFrame.from_dict({sample_id:collapsed_series}), \
             pd.DataFrame.from_dict({sample_id:feat_count_series}), \
             pd.DataFrame.from_dict({sample_id:feat_tax_series}), \
             feat_tempSeq_dict)

##main##
if __name__ == '__main__':
    if not len(sys.argv) == 5:
        print ("coverge.py countTable.csv asv_PE.fasta asv_PE.cls log")
        sys.exit()

start_time = timeit.default_timer()

#environmental variables
WD = os.environ['WD'] #work directory
RDPHOME = os.environ['RDPHOME']
RESDIR = os.environ['RESDIR']
readseqJar = os.path.join(RDPHOME, "ReadSeq.jar")
ref_db_path = os.path.join(WD, 'refset.fasta')

#Three arguments#
count_table = sys.argv[1] #countable using asv representative IDs
pe_seq_dict = read_seq(sys.argv[2]) #PE reads
read_cls = sys.argv[3] #assignment of asv PEs
log = sys.argv[4] #the name of the log file to append processing info

##load countTable
rDF = pd.read_csv(count_table, sep = ',', header=0, index_col = 0)

#Fetch tax assignments of all asv PEs prior to converging
pe_lineage_dict = get_lineages(read_cls, 0.7)

converge_start = timeit.default_timer()

#run parallel processing on seqmatch
seq_match_tmp = os.path.join(WD, 'asv_tmp')
K1_dict=run_seqmatch(seq_match_tmp, WD)

#Iterate over all samples to converge asv PEs; it's optional for parallel processing
#with multiple processors if the resources including memory are available
lineage_count_series_list = []
feature_count_series_list = []
taxonomy_series_list = []
template_mapped_seqs_dict = {} #template seqs for associated reads

###############################################################################################
#Process one sample at a time                                                                 #
#for sample_id in rDF.columns:                                                                #
#    collapsed, featCounts, featTax, templates_mapped_seqs_dict_sample = \ converge(sample_id)#
#    lineage_count_series_list.append(collapsed)                                              #
#    feature_count_series_list.append(featCounts)                                             #
#    taxonomy_series_list.append(featTax)                                                     #
#    template_mapped_seqs_dict.update(templates_mapped_seqs_dict_sample)                      #
###############################################################################################

###############################################################################
#Multi-processing by process pooling. Use caution when applying this option   #
all_completed = [] #to collect abundance for each sample when it's completed  #
with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:       #
    results = [executor.submit(converge, sample) for sample in rDF.columns]   #
    for f in concurrent.futures.as_completed(results):                        #
        all_completed.append(f.result())                                      #
for sample_rs in all_completed:                                               #
    lineage_count_series_list.append(sample_rs[0])                            #
    feature_count_series_list.append(sample_rs[1])                            #
    taxonomy_series_list.append(sample_rs[2])                                 #
    template_mapped_seqs_dict.update(sample_rs[3]) #template seqs for assoc   #
###############################################################################

# Render the abundance and taxonomy tables
lineage_abundance_table = pd.concat(lineage_count_series_list, join = 'outer',\
        axis = 1, sort=False).fillna(0).round(2)
lineage_abundance_table.to_csv(os.path.join(RESDIR, 'lineage-table.tsv'), sep='\t')
feature_abundance_table = pd.concat(feature_count_series_list, join = 'outer', \
        axis = 1, sort=False).fillna(0).round(2)
feature_abundance_table.to_csv(os.path.join(RESDIR,'feature-table.tsv'), sep='\t')
taxonomy_table = pd.concat(taxonomy_series_list, join = 'outer', axis = 1, sort=False)
taxonomy_table = taxonomy_table.fillna('-').max(axis=1, numeric_only=False)
taxonomy_table.index.name = 'feature'
taxonomy_table.to_csv(os.path.join(RESDIR, 'taxonomy-table.tsv'), sep='\t',\
                      header = ['lineage'])

# make the reference tree using the template sequences
templateIDs_all = taxonomy_table.index.to_list()
templateIDs_mapped = template_mapped_seqs_dict.keys()
templateIDs_unmapped = set(templateIDs_all).difference(set(templateIDs_mapped))
with open(os.path.join(WD, 'unmapped_templates.list'), 'w') as out:
    out.write("\n".join(templateIDs_unmapped))
    out.close()
seqid_filename = os.path.join(WD, 'unmapped_templates.list')
seq_filename = os.path.join(WD, 'templates.fasta')
fetch_refseq(RDPHOME, seqid_filename, seq_filename, ref_db_path)
#append template seqs of mapped reads
with open(os.path.join(WD, 'templates.fasta'), 'a') as out:
    for ID in templateIDs_mapped:
        out.write('>' + ID + '\n' + template_mapped_seqs_dict[ID] + '\n')

build_tree(seq_filename, WD, RESDIR)

elapsed = timeit.default_timer() - converge_start
print ('Total converge.py time: ', elapsed)
