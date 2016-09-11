"""Collection of string matching utilities for DNA sequences."""


def boyer_moore(p, p_bm, t):
    """Do Boyer-Moore matching."""
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


def naive_with_counts(p, t):
    """Return the naive match.

    Also return the number of alignments
    and the total number of character comparisons.
    """
    occurrences = []
    alignments = 0
    comparisons = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        alignments += 1
        match = True
        for j in range(len(p)):  # loop over characters
            comparisons += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, alignments, comparisons


def bm_with_counts(p, p_bm, t):
    """Do Boyer-Moore matching.

    Return the number of matches, the number of alignments
    and the total number of character comparisons.
    """
    i = 0
    occurrences = []
    alignments = 0
    comparisons = 0
    while i < len(t) - len(p) + 1:
        shift = 1
        alignments += 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            comparisons += 1
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
    return occurrences, alignments, comparisons


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
