"""Collection of string matching utilities for DNA sequences."""

from bm_preproc import BoyerMoore
from pigeonhole import Index, SubseqIndex


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


def approximate_match(p, t, n):
    """Return matches allowing n mismatches."""
    segment_length = int(len(p) / (n+1))
    all_matches = set()
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        p_bm = BoyerMoore(p[start:end], alphabet='ACGT')
        matches = boyer_moore(p[start:end], p_bm, t)
        hits = 0
        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue
            hits += 1
            mismatches = 0
            for j in range(0, start):
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)):
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            if mismatches <= n:
                all_matches.add(m-start)

    return list(all_matches), hits


def queryIndex(p, t, index, n):
    """Return matches allowing n mismatches."""
    k = index.k
    offsets = []
    for i in index.query(p):
        mismatches = 0
        for j in range(len(p[k:])):
            if p[k:][j] != t[i+k:i+len(p)][j]:
                mismatches += 1
            if mismatches > n:
                break
        if mismatches <= n:
            offsets.append(i)
    return offsets


def query_subseq(p, t, subseq_ind):
    """Return matches allowing n mismatches."""
    k = subseq_ind.k
    offsets = []
    hits = 0
    for i in subseq_ind.query(p):
        mismatches = 0
        for j in range(len(p[k:])):
            if p[k:][j] != t[i+k:i+len(p)][j]:
                mismatches += 1
            if mismatches > 2:
                break
        if mismatches <= 2:
            print("mismatches: ", mismatches)
            print("index: ", i)
            print(p[k:])
            print(t[i+k:i+len(p)])
            offsets.append(i)

    for i in range(len(subseq_ind.index)):
        hits += len(subseq_ind.query(p[i:]))
    return offsets, hits


def index_match(p, t, n):
    """Return matches allowing n mismatches."""
    segment_length = int(len(p) / (n+1))
    all_matches = set()
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        p_bm = BoyerMoore(p[start:end], alphabet='ACGT')
        matches = boyer_moore(p[start:end], p_bm, t)

        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue

            mismatches = 0
            for j in range(0, start):
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)):
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            if mismatches <= n:
                all_matches.add(m-start)

    return list(all_matches), matches
