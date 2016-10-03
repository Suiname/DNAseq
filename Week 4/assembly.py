#!/usr/bin/env python

"""Week 3 homework, algorithms for assembly.

Using DeBruijn and shortest common superstring.
"""

import itertools
from collections import defaultdict


def overlap(a, b, min_length=3):
    """Return length of longest suffix of 'a'.

    function matches a prefix of 'b' that is at least 'min_length'
    characters long.  If no such overlap exists,
    return 0.
    """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match


def scs(ss):
    """Return shortest common superstring."""
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


def all_scs(ss):
    """Return all shortest common superstrings."""
    shortest_sup = None
    all_short_sups = []
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) <= len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
            all_short_sups.append(sup)
    length = len(all_short_sups[0])
    for sup in all_short_sups:
        if len(sup) < length:
            length = len(sup)
    results = []
    for sup in all_short_sups:
        if len(sup) == length:
            results.append(sup)
    return results  # return shortest


def pick_maximal_overlap(reads, k):
    """ Return a pair of reads from the list with a
        maximal suffix/prefix overlap >= k.  Returns
        overlap length 0 if there are no such overlaps."""
    reada, readb = None, None
    best_olen = 0

    # Make index
    index = defaultdict(set)
    for read in reads:
        for i in range(len(read) - k + 1):
            index[read[i:i+k]].add(read)

    for r in reads:
        for o in index[r[-k:]]:
            if r != o:
                olen = overlap(r, o, k)
                if olen > best_olen:
                    reada, readb = r, o
                    best_olen = olen

    return reada, readb, best_olen


def greedy_scs(reads, k):
    """Return shortest common superstring using greedy algorithm."""
    read_a, read_b, olen = pick_maximal_overlap(reads, k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap(reads, k)
    return ''.join(reads)


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
