def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences


""" implement versions of the naive exact matching and Boyer-Moore algorithms that 
additionally count and return 
(a) the number of character comparisons performed and 
(b) the number of alignments tried. 
Roughly speaking, these measure how much work the two different algorithms are doing."""


def naive_with_counts(p, t):
    occurence = []
    alignments = 0
    character = 0
    for i in range(len(t)-len(p) + 1):
        alignments += 1
        match = True
        for j in range(len(p)):
            character += 1
            if not p[j] == t[i+j]:
                match = False
                break
        if match:
            occurence.append(i)
    return occurence, alignments, character


def bm_with_counts(p, p_bm, t):
    i = 0
    occurrences = []
    alignments = 0
    character = 0
    while i < len(t) - len(p) + 1:
        alignments += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            character += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, alignments, character

"""
# test 
p = 'word'
t = 'there would have been a time for such a word'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print(occurrences, num_alignments, num_character_comparisons)

p = 'needle'
t = 'needle need noodle needle'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print(occurrences, num_alignments, num_character_comparisons)

from bm_preproc import BoyerMoore
p = 'word'
t = 'there would have been a time for such a word'
lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = bm_with_counts(p, p_bm, t)
print(occurrences, num_alignments, num_character_comparisons)

p = 'needle'
t = 'needle need noodle needle'
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = bm_with_counts(p, p_bm, t)
print(occurrences, num_alignments, num_character_comparisons)
"""

# home work 

def readGenome(filename):
    """ a function that parses a DNA reference genome from a file in the FASTA format"""
    genome = ''
    with open(filename,'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome



"""
chr1 = readGenome('chr1.GRCh38.excerpt.fasta')
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, chr1)
print(num_alignments)

p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, chr1)
print(num_character_comparisons)

from bm_preproc import BoyerMoore

p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
p_bm = BoyerMoore(p)
occurrences, num_alignments, num_character_comparisons = bm_with_counts(p, p_bm, chr1)
print(num_alignments)
"""

import bisect
class Index(object):
    def __init__(self, t, k):
        ''' Create index from all substrings of size 'length' '''
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer
    
    def query(self, p):
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits



""" Write a function that, given a length-24 pattern P and given an Index object built on 8-mers, 
finds all approximate occurrences of P within T with up to 2 mismatches. Insertions and deletions are not allowed. 
Don't consider any reverse complements.
"""
def naive_2mm(p, t):
    """allows up to 2 mismatches per occurrence. 
    Unlike for the previous questions, do not consider the reverse complement here. """
    occurence = []
    for i in range(len(t)-len(p) + 1):
        match = True
        unmatch = 0
        for j in range(len(p)):
            if not p[j] == t[i+j]:
                unmatch += 1
                if unmatch > 2:
                    match = False
                    break
        if match:
            occurence.append(i)
    return occurence

# question 4 

def num_mismatch(p, t):
    mistaches = 0 
    for i in range(len(p)):
        if p[i] != t[i]:
            mistaches += 1
    return mistaches

 
def apporoximate_queryIndex(p, t, index, n):
    k = index.k
    l = len(p)
    offsets = []
    for i in range(l - k + 1):
        kmer = p[i: i+k]
        mismatches = 0 
        for j in index.query(kmer):
            if num_mismatch(p, t[j-i: j+l-i]) <= 2:
                if j-i not in offsets:
                    offsets.append(j-i)
    return offsets

# test
"""
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
chr1 = readGenome('chr1.GRCh38.excerpt.fasta')
k_mer_index = Index(chr1, 8)

r1 = apporoximate_queryIndex(p, chr1, k_mer_index, 2)
r2 = naive_2mm(p, chr1)

same = True
for r in r1:
    if r not in r2:
        same = False
print(same)

print(len(r1))
"""

# question 5

def approximate_match(p, t, n):
    from bm_preproc import BoyerMoore
    """ pigeon hole """
    segment_length = int(round(len(p)/ (n+1)))
    all_matches = set()
    num_index_hits = 0
    for i in range(n+1):
        start = i * segment_length
        end = min((i + 1) * segment_length, len(p))
        p_bm = BoyerMoore(p[start:end])
        matches = boyer_moore(p[start:end], p_bm, t)
        num_index_hits += len(matches)
        for m in matches:
            if (m < start) or (m + len(p) - start > len(t)) :
                continue

            mistaches = 0 
            for j in range(0, start):
                if p[j] != t[m-start+j]:
                    mistaches += 1
                    if mistaches > n:
                        break
            for j in range(end, len(p)):
                if p[j] != t[m+j-start]:
                    mistaches += 1
                    if mistaches > n:
                        break

            if mistaches <= n:
                all_matches.add(m-start)

    return list(all_matches), num_index_hits


p = 'GGCGCGGTGGCTCACGCCTGTAAT'
chr1 = readGenome('chr1.GRCh38.excerpt.fasta')
occurrences, num_index_hits = approximate_match(p, chr1, 2)
print(occurrences)
print(num_index_hits)


# question 6

import bisect
   
class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

"""Write a function that, given a length-24 pattern P and given a SubseqIndexobject built with k = 8 and ival = 3, 
finds all approximate occurrences of P within T with up to 2 mismatches."""

def apporoximate_querySubseqIndex(p, t, subseqindex, n):
    k = subseqindex.k
    ival = subseqindex.ival
    l = len(p)
    offsets = []
    num_index_hits = 0
    for i in range(ival):
        hits = subseqindex.query(p[i:])
        num_index_hits += len(hits)
        mismatches = 0 
        for j in hits:
            substring = t[j-i: j+l-i]
            if num_mismatch(p, substring) <= 2:
                if j-i not in offsets:
                    offsets.append(j-i)
    return offsets, num_index_hits


# test
t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
p = 'to-morrow and to-morrow '
subseq_ind = SubseqIndex(t, 8, 3)
occurrences, num_index_hits = apporoximate_querySubseqIndex(p, t, subseq_ind, 2)
print(occurrences)
print(num_index_hits)

t = open('1110.txt').read()
p = 'English measure backward'
subseq_ind = SubseqIndex(t, 8, 3)
occurrences, num_index_hits = apporoximate_querySubseqIndex(p, t, subseq_ind, 2)
print(occurrences)
print(num_index_hits)
print(t[occurrences[0]:occurrences[0]+len(p)])

# question 
chr1 = readGenome('chr1.GRCh38.excerpt.fasta')
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
subseq_ind = SubseqIndex(t = chr1, k = 8, ival = 3)
occurrences, num_index_hits = apporoximate_querySubseqIndex(p, chr1, subseq_ind, 2)

print(occurrences)
print(num_index_hits)



















