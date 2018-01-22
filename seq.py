"""
We can ask the user to supply a desired min and max fragment length.
Simulated sequencing errors and variable sequence quality can be introduced to
increase the realism of the simulation.

In a real sequencing project, fragments of the DNA to be sequenced are produced by
random processes. Thus, our program should randomly choose a substring of the input
DNA string that falls within  the specified size range. However, we need to make sure
generate enough overalpping fragments to cover the whole genome. The original shotgun
sequencing genome projects tried to achieve about eightfold coverage of the entire
original sequence: that is, each base position in the original sequence should appear in
at least eight fragments.

Next generation sequencing methods, with their shorter reads, typically work with 30-fold coverage,
while an application such as identifying rare mutations with a high degree of confidence may require
1000-fold coverage.
User input:
desired coverage value

Goal: generate random fragments from input seq
"""

import random

### User definitions
# Fragment size
fMin = 50
fMax = 100
fold = 30 # each base should appear in at least 30 fragments

def read_seq(fasta):
    seq = ''
    with open(fasta, 'r') as sequence:
        for s in sequence:
            s = s.strip()
            if len(s) > 0:
                if s[0] != '>':
                    seq += s
    return seq

seq = read_seq('homo.sapiens.mtc.fna')

def coverage(reads):
    G = len(seq)
    N = len(reads)
    L = len(''.join(reads)) / N
    cov = float(N)* (float(L)/float(G))
    return (cov)

def fragment(seq, reads):
    fSize = random.randint(fMin, fMax)
    reads.append(''.join(random.sample(seq, fSize)))
    cov = coverage(reads)

    while cov <= fold:
        fSize = random.randint(fMin, fMax)
        reads.append(''.join(random.sample(seq, fSize)))

        cov = coverage(reads)
    return reads



fragments = fragment(seq, reads = [])

def overlap(fragments):
    

contig = ''





