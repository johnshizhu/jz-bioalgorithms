import matplotlib.pyplot as plot
import gzip

'''
Args:
    seq --> a string sequence
    k   --> size of kmer
return:
    a dictionary of all kmers and frequency of kmer
'''
def kmerCounts(seq, k):
    kmerDict = {}
    for i in range(1,len(seq)-k+1):
        kmer = seq[i:i+k]
        kmerDict[kmer] = kmerDict.get(kmer,0) + 1
    return kmerDict

'''
Args:
    filename --> fasta file path
return:
    two lists, one of headers and one of sequences
'''
def loadFasta(filename):
    """ Parses a classically formatted and possibly 
        compressed FASTA file into two lists. One of 
        headers and a second list of sequences.
        The ith index of each list correspond."""
    if (filename.endswith(".gz")):
        fp = gzip.open(filename, 'r')
    else:
        fp = open(filename, 'r')
    # split at headers
    data = fp.read().split('>')
    fp.close()
    # ignore whatever appears before the 1st header
    data.pop(0)     
    headers = []
    sequences = []
    for sequence in data:
        lines = sequence.split('\n')
        headers.append(lines.pop(0))
        # add an extra "+" to make string "1-referenced"
        sequences.append('+' + ''.join(lines))
    return (headers, sequences)


# Compute a histogram of kmer-counts (i.e. how many kmers appear 1 time, 2 times, 3 times ...)
header, seq = loadFasta("data\SARS-COV-2Wuhan.fasta")
k = 6
maxcount = 50
kmers = kmerCounts(seq[0], k)
hist = [0 for i in range(maxcount)]
for kmer in kmers:
    count = kmers[kmer]
    if (count < maxcount):
        hist[count] += 1

fig = plot.figure(figsize=(10,6))
plot.title("k-mer frequency (How many k-mers, y, have a count of x)")
plot.plot([i for i in range(maxcount)], hist)
plot.show()

import random

fig = plot.figure(figsize=(10,6))
plot.title("k-mer frequency in random genomes")
for j in range(20):
    # Make a fake genome of random nucleotides
    fake = '+' + ''.join(random.choices("ACGT", k=len(seq[0])-1))
    k = 6
    maxcount = 50
    kmers = kmerCounts(fake, k)
    hist = [0 for i in range(maxcount)]
    for kmer in kmers:
        count = kmers[kmer]
        if (count < maxcount):
            hist[count] += 1
        if (count > 25):
            print(kmer, count)
    plot.plot([i for i in range(maxcount)], hist)
plot.show()

gene = {
    "ORF1a": (266, 13484),
    "ORF1ab": (266, 21556),
    "S": (21563, 25385),
    "ORF3a": (25393, 26221),
    "E": (26245, 26473),
    "M": (26523, 27192),
    "ORF6": (27202, 27388),
    "ORF7a" : (27394, 27760),
    "ORF7b": (27756, 27888),
    "ORF8": (27894, 28260),
    "N": (28274, 29534),
    "ORF10": (29558, 29675),
}

start, end = gene['S']    # Spike gene
spike = seq[0][start:end]
print(spike, len(spike))

codon = {  # RNA triplet to amino acid mapping 
    "AAA": 'K', "AAG": 'K', "AAC": 'N', "AAT": 'N',
    "AGA": 'R', "AGG": 'R', "AGC": 'S', "AGT": 'S',
    "ACA": 'T', "ACG": 'T', "ACC": 'T', "ACT": 'T',
    "ATA": 'I', "ATG": 'M', "ATC": 'I', "ATT": 'I',
    "GAA": 'E', "GAG": 'E', "GAC": 'D', "GAT": 'D',
    "GGA": 'G', "GGG": 'G', "GGC": 'G', "GGT": 'G',
    "GCA": 'A', "GCG": 'A', "GCC": 'A', "GCT": 'A',
    "GTA": 'V', "GTG": 'V', "GTC": 'V', "GTT": 'V',
    "CAA": 'Q', "CAG": 'Q', "CAC": 'H', "CAT": 'H',
    "CGA": 'R', "CGG": 'R', "CGC": 'R', "CGT": 'R',
    "CCA": 'P', "CCG": 'P', "CCC": 'P', "CCT": 'P',
    "CTA": 'L', "CTG": 'L', "CTC": 'L', "CTT": 'L',
    "TAA": '*', "TAG": '*', "TAC": 'Y', "TAT": 'Y',
    "TGA": '*', "TGG": 'W', "TGC": 'C', "TGT": 'C',
    "TCA": 'S', "TCG": 'S', "TCC": 'S', "TCT": 'S',
    "TTA": 'L', "TTG": 'L', "TTC": 'F', "TTT": 'F'
}

AminoAcid = { # Amino Acid letter Mapping
    'A': 'Alanine', 'C': 'Cysteine', 'D': 'Aspartic acid', 'E': 'Glutamic acid', 'F': 'Phenylalanine',
    'G': 'Glycine', 'H': 'Histidine', 'I': 'Isoleucine', 'K': 'Lysine', 'L': 'Leucine', 'M': 'Methionine',
    'N': 'Asparagine', 'P': 'Proline', 'Q': 'Glutamine', 'R': 'Arginine', 'S': 'Serine',
    'T': 'Theronine', 'V': 'Valine', 'W': 'Tryptophan', 'Y': 'Tyrosine', '*': 'STOP'
}

peptide = ''.join([codon[spike[i:i+3]] for i in range(0,len(spike),3)])
print(peptide)