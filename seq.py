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

# Parse fasta file
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

# Calculate coverage according to
# wikipedia.org/wiki/Coverage_(genetics)#Calculation
def coverage(reads):
    G = len(seq)
    N = len(reads)
    L = len(''.join(reads)) / N
    cov = float(N) * (float(L) / float(G))
    return (cov)

def pick_read(rloc, fSize, seq, tmpread=''):
    if (rloc+fSize) > (len(seq)-1):
        tmpread = seq[rloc:(len(seq)-1)]
        tmpread += seq[0:(rloc+fSize)-(len(seq)-1)]
        return tmpread
    else:
        tmpread = seq[rloc:(rloc+fSize)]
        return tmpread

def fragment(seq, reads):
    # Fragment size
    fSize = random.randint(fMin, fMax)
    # Random location to pick from
    rloc = random.randint(0,len(seq)-1)

    read = pick_read(rloc, fSize, seq, tmpread='')
    reads.append(read)

    cov = coverage(reads)

    while cov <= fold:
        # Fragment size
        fSize = random.randint(fMin, fMax)
        # Random location to pick from
        rloc = random.randint(0,len(seq)-1)
        # Check if out of range
        read = pick_read(rloc, fSize, seq, tmpread='')
        reads.append(read)

        cov = coverage(reads)
    return reads



fragments = fragment(seq, reads = [])
for i in range(50):
    print(fragments[i])

#def overlap(fragments):
#print(len(fragments)

#contig = ''
