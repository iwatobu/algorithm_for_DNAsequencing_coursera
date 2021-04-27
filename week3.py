def readGenome(filename):
    """ a function that parses a DNA reference genome from a file in the FASTA format"""
    genome = ''
    with open(filename,'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome



def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]

def editDistance_pt(x, y):
    # Create distance matrix , x=pattern, y=text
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first column of matrix
    # first row all 0 as we do not know when x apear in y.
    for i in range(len(x)+1):
        D[i][0] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the minimum value in the final row
    return min(D[-1])

"""
# test
p = 'GCGTATGC'  
t = 'TATTGGCTATACGGTT'
print(editDistance_pt(p, t))
"""
p = 'GCTGATCGATCGTACG'
chr1 = readGenome('chr1.GRCh38.excerpt.fasta')
print(editDistance_pt(p, chr1))

p = 'GATTTACCAGATTGAG'
chr1 = readGenome('chr1.GRCh38.excerpt.fasta')
print(editDistance_pt(p, chr1))


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists, return 0. """

    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match


def naive_overlap_map(reads, k):
	olaps = {}
	for a,b in permutations(reads, 2):
		olen = overlap(a,b, k)
		if olen > 0:
			olaps[(a,b)] = olen
	return olaps


""""
Let every k-mer in the dataset have an associated Python set object, which starts out empty.  
We use a Python dictionary to associate each k-mer with its corresponding set. 
(1) For every k-mer in a read, we add the read to the set object corresponding to that k-mer.  
If our read is GATTA and k=3, we would add GATTA to the set objects for GAT, ATT and TTA.  
We do this for every read so that, at the end, each set contains all reads containing the corresponding k-mer.  
(2) Now, for each read a, we find all overlaps involving a suffix of a.  
To do this, we take a's length-k suffix, find all reads containing that k-mer (obtained from the corresponding set) and 
call overlap(a, b, min_length=k) for each.
"""

def modified_overlap_map(reads, k):
	# k-mer dictionary
	kmer_reads = {}
	for read in reads:
		for i in range(len(read) - k + 1):
			kmer = read[i:i+k]
			if kmer not in kmer_reads.keys():
				kmer_reads[kmer] = {read}
			else:
				kmer_reads[kmer].add(read)
	
	# for each read a, find all overlaps involving a suffix of a, then do overlap 
	olaps = {}
	for a in reads:
		kmer = a[len(a)-k:len(a)]
		read_list = list(kmer_reads[kmer]) 
		for b in read_list:
			if a!=b: # no compair to itself
				olen = overlap(a,b,k)
				if olen > 0:
					olaps[(a,b)] = olen
	return list(olaps.keys())

"""
#test
reads = ['ABCDEFG', 'EFGHIJ', 'HIJABC']
print(modified_overlap_map(reads, 3))
print(modified_overlap_map(reads, 4))

reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
print(modified_overlap_map(reads, 4))
print(modified_overlap_map(reads, 5))
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



# quiz
reads, quals = readFastq('ERR.fastq')
olaps = modified_overlap_map(reads, 30)
print(len(olaps))

out = set()
for t in olaps:
	out.add(t[0])
print(len(out))









