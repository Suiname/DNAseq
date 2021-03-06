"""Run code to answer quiz questions for week 4."""

from assembly import scs, all_scs, readFastq, greedy_scs
strings = ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']
print("SCS result: ", scs(strings))
print("SCS length: ", len(scs(strings)))
all_scses = all_scs(strings)
print("SCS all result: ", all_scses)
reads, qualities = readFastq('ads1_week4_reads.fq')
genome_length = 15894
for k in range(100, 1, -1):
    full_genome = greedy_scs(reads, k)
print("Fully assembled genome", full_genome)
print("Fully assembled genome length: ", len(full_genome))
print("Actual assembled genome length: ", genome_length)
print("Number of As: ", full_genome.count("A"))
print("Number of Ts: ", full_genome.count("T"))
