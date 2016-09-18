from assembly import readGenome
from assembly import matchDistance
genome = readGenome('chr1GRCh38.excerpt.fasta')
p = 'GCGTATGC'
t = 'TATTGGCTATACGGTT'
distance, distanceMatrix = matchDistance(p, t)
print(distance)
print(distanceMatrix)
