#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

#dereplicate reference hits by their matching read id set
def derep_refset(refset):
    uniq_read_sets = []
    dereped = []
    for ref in refset:
        q_set = ref.getReadIDs()
        if not q_set in uniq_read_sets:
            dereped.append(ref)
            uniq_read_sets.append(q_set)
    return dereped

def converge_ref(refset):
    #sort refset in descending order by the number of covered bases and then by read count
    refset.sort(key=lambda x: (len(x.positions), len(x.getReadIDs())), \
                reverse=True)
    all_read_ids = set([])
    ref2read_dict = {}
    for ref in refset:
        all_read_ids = all_read_ids.union(ref.getReadIDs())
    coveredReads = set([])
    for ref in refset:
        ref_id = ref.ID
        read_ids = ref.getReadIDs()
        if not read_ids.issubset(coveredReads):
            coveredReads = coveredReads.union(ref.getReadIDs())
            ref2read_dict[ref_id] = read_ids
        if all_read_ids.issubset(coveredReads):
            break
    #ref2read_dict = remove_nonuniq_refs(ref2read_dict)
    selected = []
    for ref in refset:
        if ref.ID in ref2read_dict.keys():
            selected.append(ref)
    return selected

# recursively remove references to which there are no uniquely mapped reads
def remove_nonuniq_refs(ref2read_dict):
    read2ref_dict = {}
    for ref_id, read_ids in ref2read_dict.items():
        readIDs = ref2read_dict[ref_id]
        for read_id in read_ids:
            if not read_id in read2ref_dict:
                read2ref_dict[read_id] = set([])
            read2ref_dict[read_id].add(ref_id)
    for ref_id, read_ids in ref2read_dict.items():
        status = []
        for readID in readIDs:
            status.append(len(read2ref_dict[readID]))
        if status.count(1) == 0: #no uniquely mapped reads for this reference
            del ref2read_dict[ref_id]
            return remove_nonuniq_refs(ref2read_dict)
    return ref2read_dict

# load filtered and dereplicated blastn results into a dictionary
def get_blastn_hits(blastn, uc_filename):
    read2rep = get_id_dict(uc_filename)
    rep2asv = invert_dict(read2rep)
    hit_info = {} #to hold blast info, i.e. coordinates
    f = open(blastn, 'r')
    rc_set = set([]) #PEs in reverse-complement orientation with respect to the reference sequences (presumed to be in the positive orientation)
    while 1:
        try:
            line = f.next()
            cols = line.strip().split('\t')
            qid = cols[0]
            asv_ids = rep2asv[qid] #the actual asv_ids represented by this rep ID
            sid = cols[1].strip()
            start = int(cols[8])
            end = int(cols[9])
            asv_ids = rep2asv[qid] #the actual asv_ids represented by this rep ID
            if not sid in hit_info:
                hit_info[sid] = {}
            for asv_id in asv_ids:
                asv_id, r_end = asv_id.split('_R')#R1:0; R2:1
                r_end = int(r_end) - 1
                if not asv_id in hit_info[sid]:
                    hit_info[sid][asv_id] = [[], []]
                if start > end: # reverse-complement orientation
                    rc_set.add(asv_id)
                info = [start, end]
                info.sort()
                hit_info[sid][asv_id][r_end] = info
        except StopIteration:
            break
    return hit_info, rc_set
#generate a map of unique asvid to asvid
def get_id_dict(uc_filename):
    id_dict = {}
    f = open(uc_filename, 'r')
    while 1:
        try:
            line = f.next()
            cols = line.strip().split('\t')
            if not cols[0] == 'C':
                cols = line.strip().split('\t')
                self, rep = cols[8:10]
                if rep == '*':
                    id_dict[self] = self
                else:
                    id_dict[self] = rep
        except StopIteration:
            break
    return id_dict

def remove_blastn_subsets(Hash):
    ref_ids = Hash.keys()
    ref_IDs = ref_ids[:]
    for ref_id1 in ref_ids:
        for ref_id2 in ref_ids:
            if ref_id1 == ref_id2:
                continue
            if set(Hash[ref_id1].keys()).issubset(set(Hash[ref_id2].keys())):
                ref_IDs.remove(ref_id1)
                break
    return {ref_id:Hash[ref_id] for ref_id in ref_IDs}

def invert_dict(id_dict):#switch keys and values
    id_dict_inv = {}
    for key, value in id_dict.items():
        if not value in id_dict_inv:
            id_dict_inv[value] = []
        id_dict_inv[value].append(key)
    return id_dict_inv

# instantiate objects for template sequences
def initiate_refset(hits_info_dict):
    from Refseq import Refseq
    refset = []
    asv_list = []
    for refid in hits_info_dict.keys():
        info = hits_info_dict[refid]
        asv_ids = info.keys()
        asv_list.extend(asv_ids)
        if not refid in refset:
            myref = Refseq(refid)
        for asv_id in asv_ids:
            myref.addPEread(asv_id, info[asv_id][0], info[asv_id][1])
        refset.append(myref)
    asv_list = list(set(asv_list))
    return refset, asv_list
