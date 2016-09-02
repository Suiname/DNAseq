"""Collection of string matching utilities for DNA sequences."""


def naive(p, t):
    """Return the naive match of a substring to a string."""
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences


def reverseComplement(s):
    """Return the reverse complement of a string."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


def readGenome(filename):
    """Return a sequence from a fasta file."""
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def readFastq(filename):
    """Return a tuple of sequences and qualities from a fastq file."""
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip()  # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


# implement a version of naive matching that is strand aware
def sa_naive(p, t):
    """Return a naive string match that considers the reverse complement."""
    occurrences = []
    if reverseComplement(p) == p:  # if the reverse complement is the same,
        return naive(p, t)  # just return the naive match
    occurrences = naive(p, t)  # otherwise, set result to be the naive match
    occurrences += naive(reverseComplement(p), t)  # then add in the complement
    return occurrences


def naive_2mm(p, t):
    """Return a naive string match allowing up to 2 mismatches."""
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        mismatch = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mismatch += 1  # increment mismatch
                if mismatch > 2:  # mismatches > 2 makes it a mismatch
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences
