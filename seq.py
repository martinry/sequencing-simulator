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

seq = read_seq('homo10.fna')

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

def find_kmers(reads, k):
    kmers = []
    for read in fragments:
        s = [read[i:i+k] for i in range(len(read)) if len(read[i:i+k]) == k]
        kmers.extend(s)
    kmers = set(kmers)
    return list(kmers) # return unique kmers

# k: kmer size
kmers = find_kmers(fragments, 4)

de_bruijn = {x:[] for x in kmers}

for m in kmers:
    for k in kmers:
        #print(kmers[0], k)
        #print(kmers[0][1:], k[0:-1])
        if m[1:] == k[0:-1]:
            de_bruijn[m].append(k)


for k,v in de_bruijn.items():
    #v = [(x[0], x[1:]) for x in v]
    print('%s(%s):\t%s' % (k[0],k[1:],v))
#print(len(kmers))









#contig = ''
