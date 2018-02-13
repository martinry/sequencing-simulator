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

seq = read_seq('h10.fna')

#seq = 'AATGCCGTACGTAGGGTAATATATGACCA'

# Calculate coverage according to
# wikipedia.org/wiki/Coverage_(genetics)#Calculation
def coverage(reads):
    G = len(seq)    # Length of genome sequence
    N = len(reads)  # Number of reads
    L = len(''.join(reads)) / N # Average read length
    cov = float(N) * (float(L) / float(G))
    return (cov)

def pick_read(rloc, fSize, seq, tmpread=''):
    # Account for circular genome:
    # if the 'picked' length is greater than the length of the genome seq
    # We have to continue picking in the 'beginning' of the seq
    if (rloc+fSize) > (len(seq)):
        tmpread = seq[rloc:(len(seq))]
        tmpread += seq[0:(rloc+fSize)-(len(seq))]
        return tmpread
    
    else:
        tmpread = seq[rloc:(rloc+fSize)]
        return tmpread

def fragment(seq, reads):
    # Fragment size
    fSize = random.randint(fMin, fMax)
    
    # Random position to pick from
    rloc = random.randint(0,len(seq)-1)

    # Pick a read of length fSize at position rloc
    read = pick_read(rloc, fSize, seq, tmpread='')
    reads.append(read)

    # Calculate current coverage
    cov = coverage(reads)

    # Repeat until coverage requirement has been met
    while cov <= fold:
        # Fragment size
        fSize = random.randint(fMin, fMax)
        
        # Random location to pick from
        rloc = random.randint(0,len(seq)-1)
        
        # Pick a read of length fSize at position rloc
        read = pick_read(rloc, fSize, seq, tmpread='')
        reads.append(read)

        cov = coverage(reads)
    return set(reads)

print('Generating fragments...')
fragments = fragment(seq, reads = [])

def find_kmers(reads, k):
    kmers = []
    for read in fragments:
        s = [read[i:i+k] for i in range(len(read)) if len(read[i:i+k]) == k]
        kmers.extend(s)
    kmers = set(kmers)
    return list(kmers) # return unique kmers

# the lower the kmer size, the more kmers in total
# but fewer unique kmers

print('Generating kmers...')
# k: kmer size
kmers = find_kmers(fragments, 7)


graph = {}

print('Generating graph...')
for k in kmers:
    from_node = k[1:]
    to_nodes = [m for m in kmers if from_node in m[0:-1]]
    graph[k] = to_nodes
    
import networkx as nx
import matplotlib.pyplot as plt
#%matplotlib inline


# Plot directed graph, use max kmer 8, max genome 10
def draw_nw(g):
    plt.figure(figsize=(20,20), dpi=50)
    nx.draw_networkx(g, arrows=True, with_labels=False, node_size=200, alpha=0.5, font_size=20)

print('Plotting...')
g = nx.DiGraph(graph)
#draw_nw(g)

def draw_de_bruijn_graph(g):
    plt.figure(figsize=(20,20), dpi=80)
    nx.draw_networkx(
        g, pos=nx.circular_layout(g),
        node_shape='o', alpha=0.5, with_labels=True, node_size=200, font_size=15,
        edge_color='#555555', width=0.5
    )
    nx.draw_networkx_edge_labels(
        g, pos=nx.circular_layout(g), 
        edge_labels=nx.get_edge_attributes(g, 'weight'),
        font_size=24, label_pos=0.25, rotate=False
    )
    plt.axis('off')
    plt.show()



#g = nx.DiGraph(graph)
#draw_de_bruijn_graph(g)

max_path_length = 0
longest_path_node = ''

for n in g.nodes():
    G = nx.bfs_tree(g, n)
    lpl = nx.dag_longest_path_length(G)

    if nx.dag_longest_path_length(G) > max_path_length:
        max_path_length = lpl
        longest_path_node = n

G = nx.bfs_tree(g, longest_path_node)
draw_nw(G)

print(g.size())
print(g.order())



"""
g = nx.DiGraph(graph)



nx.draw(g)
plt.savefig("simple_path.png") # save as png
plt.show() # display

try:
  path = nx.dag_longest_path(G)
  print(path)
  # ['a', 'b', 'c', 'd']

  print(len(path) - 1)
  # 3
except nx.exception.NetworkXUnfeasible: # There's a loop!
  print("The graph has a cycle")
"""
