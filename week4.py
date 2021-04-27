def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match



import itertools

def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest



# modify, return the number of different scs 
def scs_list(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    scs_ss = set()
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) <= len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
            scs_ss.add(shortest_sup)

    # find the shortest length
    scs_ss = list(scs_ss)
    l = len(scs_ss[0])
    for s in scs_ss:
        if len(s) < l:
            l = len(s)
    
    result = set()
    for s in scs_ss:
        if len(s) == l:
            result.add(s) 

    return result 

    """
    for s in scs_ss:
        if len(s) > l:
            scs_ss.remove(s)
    return scs_ss"""
    # why it doesn't work? 
    # You are not permitted to remove elements from the list while iterating over it using a for loop.


"""
# test examples 
strings = ['ABC', 'BCA', 'CAB']
print(scs(strings))
print(scs_list(strings))

strings = ['GAT', 'TAG', 'TCG', 'TGC', 'AAT', 'ATA']
print(scs(strings))
print(scs_list(strings))


#quiz
ss = ['CCT','CTT','TGC','TGG','GAT','ATT']
print(scs(ss))
print(len(scs(ss)))
print(scs_list(ss))
print(len(scs_list(ss)))
"""


def readFastq(filename):
    """ a function that parses the read and quality strings from a FASTQ file containing sequencing reads"""
    sequences = []
    qualities = []
    with open(filename, 'r') as f:
        while True: 
            f.readline() # skip name line
            seq = f.readline().rstrip()
            f.readline() # skip place holder line 
            q = f.readline().rstrip()
            if len(seq) ==0:
                break 
            sequences.append(seq)
            qualities.append(q)
    return sequences, qualities



def pick_maximal_overlap(reads,k):
    reada, readb = None, None
    best_olen = 0
    for a, b in itertools.permutations(reads, 2):
        olen = overlap(a, b, min_length = k)
        if olen > best_olen:
            reada, readb = a, b
            best_olen = olen
    return reada, readb, best_olen


def greedy_scs(reads, k):
    reada, readb, olen = pick_maximal_overlap(reads, k)
    while olen > 0:
        reads.remove(reada)
        reads.remove(readb)
        reads.append(reada + readb[olen:])
        reada, readb, olen = pick_maximal_overlap(reads, k)
    return "".join(reads)




""" 
Assemble these reads using one of the approaches discussed, such as greedy shortest common superstring.  
Since there are many reads, you might consider ways to make the algorithm faster, 
such as the one discussed in the programming assignment in the previous module.
"""

def find_kmer_readset(reads, k):
    """ find all kmer:readset dictionary
    """
    kmer_reads = {}
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            if kmer not in kmer_reads.keys():
                kmer_reads[kmer] = {read}
            else:
                kmer_reads[kmer].add(read)
    return kmer_reads


def find_max_kmer_overlap(kmer_reads,k):
    """ take the output from find_kmer_readset(), and return the max overlap reada, readb, and olen"""

    # max overlap in each kmer readset
    olaps = {}
    for reads in kmer_reads.values():
        reada, readb, olen = pick_maximal_overlap(list(reads), k)
        olaps[(reada, readb)] = olen
        

    # max overlap accross the kmer set, the first item
    readpair = max(olaps, key=olaps.get)
    return readpair[0], readpair[1], olaps[readpair]

    


def fast_greedy_scs(reads, k):
    """
    get the kmer set, 
    get the max overlap
    loop 
    """
    kmer_reads = find_kmer_readset(reads, k)
    read_a, read_b, olen = find_max_kmer_overlap(kmer_reads, k)
    while olen > 0:
        print(olen)
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        kmer_reads = find_kmer_readset(reads, k)
        read_a, read_b, olen = find_max_kmer_overlap(kmer_reads, k)
    """
    k -= 1
    while k > 0:
        print(k)
        kmer_reads = find_kmer_readset(reads, k)
        read_a, read_b, olen = find_max_kmer_overlap( kmer_reads, k)
        while olen > 0:
            reads.remove(read_a)
            reads.remove(read_b)
            reads.append(read_a + read_b[olen:])
            kmer_reads = find_kmer_readset(reads, k)
            read_a, read_b, olen = find_max_kmer_overlap(kmer_reads, k)
        k -= 1
    
    return "".join(reads)
    """
    return "".join(reads)
    #return scs(reads)
    
# 15,894  base 

# quiz
ss, _ = readFastq('ads1_week4_reads.fq')
print(len(ss))
# print(ss[:3])
g1 = fast_greedy_scs(ss,30)
print(len(g1))
print(g1.count('A'))
print(g1.count('T'))











