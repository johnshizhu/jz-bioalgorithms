# Genome Assembly
Uinsg Graphical approach to assemble a genome from fragments. 

### Kmers Unique
```python
def kmersUnique(seq, k):
    """ extracts all *k*-mers from *seq* string, while appending a
        unique subscript to each repeated k-mer """
    kmers = sorted([seq[i:i+k] for i in range(len(seq)-k+1)])  # trick is to sort them first making repeats adjacent
    for i in range(1,len(kmers)):
        if (kmers[i] == kmers[i-1][0:k]):           # check adjacent k-mers
            t = kmers[i-1].find('_')
            if (t >= 0):                            # more than 2 repeats
                n = int(kmers[i-1][t+1:]) + 1
                kmers[i] = kmers[i] + "_" + str(n)
            else:                                   # first repeat
                kmers[i-1] = kmers[i-1] + "_1"
                kmers[i] = kmers[i] + "_2"
    return kmers
```
