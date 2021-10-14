#!/usr/bin/env python
#return only doubleton and multipton clusters
#Implement centroid clustering on sorted sam pairwise alignment
#query ids are used to generate centroids, which 'recruit' hits within the 
#distance radius cutoff: defaul 0.03

import sys
from cigar import Cigar

def get_pct_id(cigar_tuple):
    m = 0 #matches
    total = 0
    for seg in cigar_tuple:
        count, Type = seg
        total += count
        if Type == '=':
            m += int(count)
    try:
        identity = m/total
    except ZeroDivisionError:
        identity = 0
    return identity

def get_dist_dict(f): #generate clusters by centroid approach
    clusters = {} #dict of clusters
    covered = set([]) #previously clustered hit IDs
    centroids = set([]) #set of centroid IDs
    order = [] #track the cluster order
    lines = open(f, 'r').readlines()
    for line in lines:
        if line[0] == '@':
            continue
        cols = line.strip().split()
        id1 = cols[0].strip()
        id2 = cols[2].strip()
        if id1 == id2:
            continue
        else:
            cigarStr = cols[5]
            c = Cigar(cigarStr)
            C = list(c.items())
            dist = 1 - get_pct_id(C)
        if dist >= cutoff:
            continue
        else:
            qid = id1
            sid = id2
            if (sid in covered) or (sid in centroids):
                continue
            if qid in covered:
                continue
            if not qid in centroids:
                clusters[qid] = [qid, sid]
                centroids.add(qid)
                order.append(qid)
            else:
                clusters[qid].append(sid)
        covered.add(sid)
    print ('centroid', len(centroids))
    print (len(covered.union(centroids)), 'in clusters')
    return clusters, order

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print ('parse_minimap2.py aligned.sam outname cutoff[0-1]')
        sys.exit()
    f = sys.argv[1]#alignment in sam format
    outname = sys.argv[2]
    cutoff = float(sys.argv[3]) #0-1 distance cutoff
    print (f)
    count = 0
    exit_set = [] #checking
    out = open(outname, 'w')
    clusters, order = get_dist_dict(f)
    for i, cID in enumerate(order):
        out.write(str(i) + '\t' + ','.join(clusters[cID]) + '\n')
    out.close()
