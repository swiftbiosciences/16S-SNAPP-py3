#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

## to parse the classifier results into a hash of counts keyed by lineage name
def get_lineages(tax_filename, CONF):
    import string
    Hash = {}
    lines = open(tax_filename.strip(), 'r').readlines()
    for row in lines:
        cols = row.strip().split('\t')
        ID = cols[0]
        ID = ID.split(';')[0]
        levels = cols[2:]
        i = 0
        lineage = []
        while i < len(levels):
            level = levels[i:i+3]
            name, rank, conf = level
            name = rank[0] + '__' + name.replace('"', '')
            if float(conf) < CONF:
                break
            else:
                lineage.append(name.replace('"', ''))
            i += 3
        lineage = string.join(lineage, ';').strip()
        if lineage == '':#rare cases the sequence can't be classified to domain level
            lineage = 'Unclassified'
        Hash[ID] = lineage
    return Hash

def get_best_tax(ref_tax, read_taxs):#choose the better between reftax and readtax
    taxs = {}
    for tax in read_taxs:
        level = tax.count(';')
        if not level in taxs:
            taxs[level] = []
        taxs[level].append(tax)
    resolution = taxs.keys()
    resolution.sort()
    resolution.reverse()
    top_lineage= taxs[resolution[0]][0]
    if ref_tax.count(';') >= resolution[0]:
        top_lineage = ref_tax
    return  top_lineage  #the most resolved lineage

def add_lineages(sample_id, tax_dict, refset):#refset pass through to add lineage information
    ref_tax_file = sample_id + '.cls'
    ref_tax_hash = get_lineages(ref_tax_file, 0.7)
    for ref_id in refset.keys():
        read_ids = list(refset[ref_id].getReadIDs())
        ref_lineage = ref_tax_hash[ref_id]
        read_lineage_list = list(set([tax_dict[ID][0] for ID in read_ids]))
        best_lineage = get_best_tax(ref_lineage, read_lineage_list)
        refset[ref_id].addAssign(best_lineage)
    return refset
