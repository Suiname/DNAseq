from assembly import readGenome
from assembly import matchDistance
genome = readGenome('chr1GRCh38.excerpt.fasta')
p = 'GCTGATCGATCGTACG'
distance = matchDistance(p, genome)
print("Distance Match of 𝙶𝙲𝚃𝙶𝙰𝚃𝙲𝙶𝙰𝚃𝙲𝙶𝚃𝙰𝙲𝙶 is: ", distance)
newp = 'GATTTACCAGATTGAG'
distance = matchDistance(newp, genome)
print("Distance Match of 𝙶𝙰𝚃𝚃𝚃𝙰𝙲𝙲𝙰𝙶𝙰𝚃𝚃𝙶𝙰𝙶 is: ", distance)
