#!/usr/bin/env python
#
import sys
import os
import pandas as pd
import json

def get_idmap(f):#replace sequence ids with integers
    Map = {}
    counter = 0
    for line in open(f, 'r').readlines():
        if line[0] == '>':
            ID = line.strip().replace('>', '')
            counter += 1
            newid = 'C' + (6-len(str(counter)))*'0' + str(counter)
            Map[ID] = newid
    return Map

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print ("get_OTU_table.py wd all_cons.fasta all_cons.clr outname")
        sys.exit()

wd = sys.argv[1] #work directory
fasta = sys.argv[2] #abundance sorted concensus fasta file
clr = sys.argv[3] #cluster file in RDPTools format
outname = sys.argv[4] #name prefix
RESDIR = os.environ['RESDIR']
outname = os.path.join(RESDIR, outname)

count_table = {}
idmap = get_idmap(fasta)

#save idmap into a json file in case needed in trouble shooting
with open("all_cons_idmap.json", "w") as outfile:
    json.dump(idmap, outfile)

#open count old table dictionary
with open('all_counts.json') as json_file:
    count_table  = json.load(json_file)

#load taxonomy table dictionary
with open('all_assignments.json') as json_file:
    tax_dict  = json.load(json_file)

#convert the cluster file into the cluster dictionar
cluster_map = {}
for line in open('all_cons.clr', 'r').readlines():
    mems = line.strip().split('\t')[1].split(',')
    cluster_map[mems[0]] = mems

#convert consensus ids to their centroid id of their respective clusters
to_centroid = {}
for cluster in list(cluster_map.values()):
    for ID in cluster:
        to_centroid[ID] = cluster[0]

#make a new abundance table with simplified centroid ids
new_tax = {}
new_dict = {}
for sample_name in list(count_table.keys()):
    new_dict[sample_name] = {}
    IDs = list(count_table[sample_name].keys())
    for ID in IDs:
        count = int(count_table[sample_name][ID])
        if ID.find('asv') == 0:#asv IDs strip sample and size info
            cID = ID.split(';')[0]
            tax = tax_dict[cID]
        else:
            try:
                cID = to_centroid[ID]
            except KeyError:
                cID = ID
            tax = tax_dict[cID]
            cID = idmap[cID] #simplified ID
        if not cID in new_dict[sample_name].keys():
            new_dict[sample_name][cID] = count
        else:
            new_dict[sample_name][cID] += count
        if not cID in new_tax.keys():
            new_tax[cID] = tax

#output the abundance table in simplified IDs
df = pd.DataFrame(new_dict).fillna(0)
df.to_csv(outname + '_count.txt', sep = '\t')

tax_out = open(outname + '_taxonomy.txt', 'w')
for f, t in new_tax.items():
    tax_out.write(f + '\t' + t + '\n')
