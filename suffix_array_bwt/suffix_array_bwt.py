def argsort(input):
    return sorted(range(len(input)), key=input.__getitem__)

def suffixArray(string):
    return sorted(range(len(string)), key=lambda x: string[x:])

def findFirst(pattern, text, suffixarray):
    lo, hi = 0, len(text)
    while (lo < hi):
        middle = (lo+hi)//2
        if text[suffixarray[middle]:] < pattern:
            lo = middle + 1
        else:
            hi = middle
    return lo

def findLast(pattern, text, suffixarray):
    lo, hi = 0, len(text)
    while (lo < hi):
        middle = (lo+hi)//2
        if text[suffixarray[middle]:suffixarray[middle]+len(pattern)] <= pattern:
            lo = middle + 1
        else:
            hi = middle
    return lo

def computeLCP(text, suffixarray):
    '''
    args
        text --> text to find
        suffixarray --> finding the longest common prefix from this
    return
        longest comon prefix
    '''
    m = len(text)
    lcp = [0 for i in range(m)]
    for i in range(1,m):
        u = suffixarray[i-1]
        v = suffixarray[i]
        n = 0
        while text[u] == text[v]:
            n += 1
            u += 1
            v += 1
            if (u >= m) or (v >= m):
                break
        lcp[i] = n
    return lcp

def BWT(t):
    '''
        Create btw of t
    '''
    # create a sorted list of all cyclic suffixes of t
    rotation = sorted([t[i:]+t[:i] for i in range(len(t))])
    # concatenate the last symbols from each suffix
    return ''.join(r[-1] for r in rotation)

def BWTfromSuffixArray(text, suffixarray):
    return ''.join(text[i-1] for i in suffixarray)

def inverseBWT(bwt):
    # initialize the table from t
    table = ['' for c in bwt]
    for j in range(len(bwt)):
        #insert the BWT as the first column
        table = sorted([c+table[i] for i, c in enumerate(bwt)])
    #return the row that ends with ‘$’

def FMIndex(bwt):
    fm = [{c: 0 for c in bwt}]  # a list of dictionaries
    for c in bwt:
        row = {symbol: count + 1 if (symbol == c) else count for symbol, count in fm[-1].items()}
        fm.append(row)
    offset = {}
    N = 0
    for symbol in sorted(row.keys()):
        offset[symbol] = N
        N += row[symbol]
    return fm, offset

def recoverSuffix(i, BWT, FMIndex, Offset):
    suffix = ''
    c = BWT[i]
    predec = Offset[c] + FMIndex[i][c]
    suffix = c + suffix
    while (predec != i):
        c = BWT[predec]
        predec = Offset[c] + FMIndex[predec][c]
        suffix = c + suffix
    return suffix

def findBWT(pattern, FMIndex, Offset):
    lo = 0
    hi = len(FMIndex) - 1
    for symbol in reversed(pattern):
        lo = Offset[symbol] + FMIndex[lo][symbol]
        hi = Offset[symbol] + FMIndex[hi][symbol]
    return lo, hi

def mergeBWT(bwt1, bwt2):
    interleave = [(c, 0) for c in bwt1] + [(c, 1) for c in bwt2]
    passes = min(len(bwt1), len(bwt2))
    for p in range(passes):
        i, j = 0, 0
        nextInterleave = []
        for c, k in sorted(interleave, key=lambda x: x[0]):
            if (k == 0):
                b = bwt1[i]
                i += 1
            else:
                b = bwt2[j]
                j += 1
            nextInterleave.append((b, k))
        if (nextInterleave == interleave):
            break
        interleave = nextInterleave
    return ''.join([c for c, k in interleave])