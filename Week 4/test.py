"""Run code to answer quiz questions for week 4."""

from assembly import scs, all_scs, readFastq, dykSuperstring
strings = ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']
print("SCS result: ", scs(strings))
all_scses = all_scs(strings)
print("SCS all result: ", all_scses)
# for ss in all_scses:
#     print("Length of SS ", ss, ": ", len(ss))
# strings = ['ABC', 'BCA', 'CAB']
# print(all_scs(strings))

reads, qualities = readFastq('ads1_week4_reads.fq')
genome_length = 15894
read_length = len("GTCCAGCAGAGCAAGTGATGCGAGAGCTGCCCATCCTCCAACCAGCATGCCCCTAGACATT\
GACACTGCATCGGAGTCAGGCCAAGATCCGCAGGACAGT")
full_genome = dykSuperstring(reads)
print("Fully assembled genome", full_genome)
print("Fully assembled genome length: ", len(full_genome))
print("Actual assembled genome length: ", genome_length)
print("Number of As: ", full_genome.count("A"))
print("Number of Ts: ", full_genome.count("T"))
