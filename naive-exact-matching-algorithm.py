import random


# This first version of the matching algorithm is the slowest and most iterative


def read_genome(filename):
    genome = ''
    with open(filename, 'r') as file:
        for line in file:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


# Parsing a FASTQ file
# 1st line: name of read ,optional: experiment name,machine used,cluster location..
# 2nd line: base sequence called by caller
# 3rd line: ignored
# 4th line: encoded base quality for each base

def read_fastq(filename):
    sequences = []
    qualities = []
    with open(filename) as input_file:
        while True:
            input_file.readline()
            seq = input_file.readline().rstrip()
            input_file.readline()
            qual = input_file.readline().rstrip()
            if len(seq) == 0:
                break
            else:
                sequences.append(seq)
                qualities.append(qual)
    return sequences, qualities


# When the DNA is devided into fragments ,the orientation of the read is not known.
# We need to check the reverse complement to ensure if two strands match

def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse = ''
    for base in dna_sequence:
        reverse = complement[base] + reverse
    return reverse


# Naive exact matching algorithm


def naive(p, t):  # Looking for pattern p in string t
    occurences = []
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):  # For each position in t ,check if p appears
            if p[j] != t[i + 1]:
                match = False
                break
        if match:
            occurences.append(i)
    return occurences


# Modification of the naive algorithm so it is strand aware:

def naive_rc(p, t):
    occurences = []
    rc = reverse_complement(p)
    for i in range(len(t) - len(p) + 1):
        if p == rc:
            match = True
            for j in range(len(p)):
                if p[j] != t[i + 1]:
                    match = False
                    break
            if match:
                occurences.append(i)
        else:
            match = True
            match_reverse = True
            for j in range(len(p)):
                if p[j] != t[i + j]:
                    match = False
                    break
            for j in range(len(rc)):
                if rc[j] != t[i + j]:
                    match_reverse = False
                    break
            if match or match_reverse:
                occurences.append(i)
    return occurences


def generate_reads(genome, num_reads, read_len):
    reads = []
    for _ in range(num_reads):
        start = random.randint(0, len(genome) - read_len)
        reads.append(genome[start: start + read_len])
    return reads


genome1 = read_genome('lambda_virus.fa')
reads1 = generate_reads(genome1, 100, 100)

# Test 1
# Testing naive exact matching algorithm on lambda virus genome and manually(fake) generating the reads

num_match = 0
for read in reads1:
    matches = naive(read, genome1)
    if len(matches) > 0:
        num_match += 1
print('%d / %d reads matched exactly!' % (num_match, len(reads1)))

# Test 2
# testing on phix.fa reads ,the real reads,not fake ones
# phix_reads, _ = read_fastq('phix.fa')

phix_reads, _ = read_fastq('C:\\Users\\Dusan\\PycharmProjects\\read-alignment-problem\\ERR266411_1.first1000.fastq')
num_matched = 0
for r in phix_reads:
    matches = naive(r, genome1)
    if len(matches) > 0:
        num_matched += 1
print('%d / %d reads matched the genome' % (num_matched, len(phix_reads)))

# Test for strand aware
num_matched = 0
for read in phix_reads:
    read = read[:30]  # only first 30 bases since this is a real long read
    match1 = naive_rc(read, genome1)
    if len(match1) > 0:
        num_matched += 1
print('%d / %d reads matched the genome exactly!' % (num_matched, len(phix_reads)))
