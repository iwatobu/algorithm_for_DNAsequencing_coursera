def naive(p, t):
	""" naive exact matching; pattern and text """
	occurence = []
	for i in range(len(t)-len(p) + 1):
		match = True
		for j in range(len(p)):
			if not p[j] == t[i+j]:
				match = False
				break
		if match:
			occurence.append(i)
	return occurence

def reverseComplement(s):
	""" function takes a DNA string and returns its reverse complement"""
	complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
	t = ''
	for base in s:
		t = complement[base] + t
	return t

def readGenome(filename):
	""" a function that parses a DNA reference genome from a file in the FASTA format"""
	genome = ''
	with open(filename,'r') as f:
		for line in f:
			if not line[0] == '>':
				genome += line.rstrip()
	return genome

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


# homework 
def naive_with_rc(p, t):
	""" naive exact matching algorithm that is strand-aware. 
	That is, instead of looking only for occurrences of P in T, additionally look for occurrences of thereverse complement of P in T. 
	If P is ACT, your function should find occurrences of both ACTand its reverse complement AGT in T"""
	rc = reverseComplement(p) 
	if p == rc:
		oc = naive(p, t)
	else:
		oc = naive(p, t) + naive(rc, t)
	return oc


# test example 1
p = 'CCC'
ten_as = 'AAAAAAAAAA'
t = ten_as + 'CCC' + ten_as + 'GGG' + ten_as
occurrences = naive_with_rc(p, t)
print(occurrences)


# test example 2
p = 'CGCG'
t = ten_as + 'CGCG' + ten_as + 'CGCG' + ten_as
occurrences = naive_with_rc(p, t)
print(occurrences)

# test example 3

phix_genome = readGenome('phix.fa')
occurrences = naive_with_rc('ATTA', phix_genome)
print('offset of leftmost occurrence: %d' % min(occurrences))
print('# occurrences: %d' % len(occurrences))


## home work 
lambda_virus_genome = readGenome('lambda_virus.fa')
p = 'AGGT'
oc = naive_with_rc(p, lambda_virus_genome)
print(len(oc))

p = 'TTAA'
oc = naive_with_rc(p, lambda_virus_genome)
print(len(oc))

p = 'ACTAAGT'
oc = naive_with_rc(p, lambda_virus_genome)
print(min(oc))

p = 'AGTCGA'
oc = naive_with_rc(p, lambda_virus_genome)
print(min(oc))


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


# test new function 
print(naive_2mm('ACTTTA', 'ACTTACTTGATAAAGT'))

p = 'CTGT'
ten_as = 'AAAAAAAAAA'
t = ten_as + 'CTGT' + ten_as + 'CTTT' + ten_as + 'CGGG' + ten_as
occurrences = naive_2mm(p, t)
print(occurrences)

phix_genome = readGenome('phix.fa')
occurrences = naive_2mm('GATTACA', phix_genome)
print('offset of leftmost occurrence: %d' % min(occurrences))
print('# occurrences: %d' % len(occurrences))

# 

print(len(naive_2mm('TTCAAGCC', lambda_virus_genome)))

p = 'AGGAGGT'
oc = naive_2mm(p, lambda_virus_genome)
print(min(oc))


def phred33ToQ(qual):
	return ord(qual) - 33

seqs, quals = readFastq('ERR.fastq')

#print(seqs[:5])
#print(quals[:5])
#print(len(seqs[0])) #100
#print(quals[0])
# build a quality matirx 
m = [ [0 for i in range(100)] for j in range(len(quals))]

#print(m[2][1:5])
for i in range(100): # column 
	for j in range(len(quals)): # row 
		qual = quals[j][i]
		m[j][i] = phred33ToQ(qual)



import numpy as np
m = np.matrix(m)
mean_by_position = m.mean(0)
print(mean_by_position)
print(mean_by_position.min())

# get the index
print(np.argmin(mean_by_position))


