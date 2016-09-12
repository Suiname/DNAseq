"""Testing."""

from pigeonhole import SubseqIndex
from better_matching import readGenome
from better_matching import query_subseq
from better_matching import approximate_match
from better_matching import naive_2mm

t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
p = 'to-morrow and to-morrow '
subseq_ind = SubseqIndex(t, 8, 3)
occurrences, num_index_hits = query_subseq(p, t, subseq_ind)
print("Index: ", subseq_ind.index)
print("Occurrences: ", occurrences)
print("num_index_hits: ", num_index_hits)

t = open('1110.txt.utf-8').read()
print("Length of T: ", len(t))
p = 'English measure backward'
subseq_ind = SubseqIndex(t, 8, 3)
occurrences, num_index_hits = query_subseq(p, t, subseq_ind)
print("Occurrences: ", occurrences)
print("num_index_hits: ", num_index_hits)
p = 'GGCGCGGTGGCTCACGCCTGTAAT'

genome = readGenome('chr1.GRCh38.excerpt.fasta')
subseq_ind = SubseqIndex(genome, 8, 3)
occurrences, hits = query_subseq(p, genome, subseq_ind)
print("Occurrences: ", sorted(occurrences))
print("num_index_hits: ", hits)

occurrences, hits = approximate_match(p, genome, 2)
print("Occurrences: ", sorted(occurrences))
print("num_index_hits: ", hits)

occurrences = naive_2mm(p, genome)
print("Occurrences: ", sorted(occurrences))
print(genome[84775:84775+24])
print(p)
