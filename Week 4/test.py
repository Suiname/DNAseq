"""Run code to answer quiz questions for week 4."""

from assembly import scs, all_scs, readFastq, greedy_scs
strings = ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']
print("SCS result: ", scs(strings))
print("SCS all result: ", all_scs(strings))
reads, qualities = readFastq('ads1_week4_reads.fq')
genome_length = 15894
full_genome = greedy_scs(reads, genome_length)
print("Fully assembled genome", full_genome)
print("Fully assembled genome length: ", len(full_genome))
print("Number of As: ", full_genome.count("A"))
print("Number of Ts: ", full_genome.count("T"))
