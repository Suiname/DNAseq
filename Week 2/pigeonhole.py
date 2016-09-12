#!/usr/bin/env python
"""kmer_index.py: A k-mer index for indexing a text."""

import bisect


__author__ = "Jason Tham"


class Index(object):
    """ Holds a substring index for a text T """

    def __init__(self, t, k):
        """ Create index from all substrings of t of length k """
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        """ Return index hits for first k-mer of p """
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

class ApproxIndex(object):
    """Holds a substring index for a text T.

    This class authored by Ben Langmead
    """

    def __init__(self, t, k):
        """Create index from all substrings of t of length k."""
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        """Return index hits for first k-mer of p."""
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            mismatch = 0
            for j in range(len(self.index[i][0])):
                if self.index[i][0][j] != kmer[j]:
                    mismatch += 1
                if mismatch > 2:
                    break
            if mismatch <= 2:
                print("self.index[0]:", self.index[i][0])
                print("Kmer:", kmer)
                print("self.index[1]:", self.index[i][1])
                hits.append(self.index[i][1])
            i += 1
        return hits


class SubseqIndex(object):
    """Holds a subsequence index for a text T."""

    def __init__(self, t, k, ival):
        """Create index from all subsequences.

        Subsequences consist of k characters spaced ival positions apart.
        E.g., SubseqIndex("ATAT", 2, 2) extracts ("AA", 0) and ("TT", 1).
        """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            # add (subseq, offset)
            self.index.append((t[i:i+self.span:ival], i))
        self.index.sort()  # alphabetize by subseq

    def query(self, p):
        """Return index hits for first subseq of p."""
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
