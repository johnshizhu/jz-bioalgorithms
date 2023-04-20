import numpy as np

def findLCS(v, w):
    #finding the LCS
    score = np.zeros((len(v)+1,len(w)+1), dtype="int32")
    backt = np.zeros((len(v)+1,len(w)+1), dtype="int32")
    for i in range(1,len(v)+1):
        for j in range(1,len(w)+1):
            # find best score at each vertex
            if (v[i-1] == w[j-1]):  # test for a match ("diagonal street")
                score[i,j], backt[i,j] = max((score[i-1,j-1]+1,3), (score[i-1,j],1), (score[i,j-1],2))
            else:
                score[i,j], backt[i,j] = max((score[i-1,j],1), (score[i,j-1],2))
    return score, backt

def LCS(b,v,i,j):
    # Backtracking function
    if ((i == 0) and (j == 0)):
        return ''
    elif (b[i,j] == 3):
        return LCS(b,v,i-1,j-1) + v[i-1]
    elif (b[i,j] == 2):
        return LCS(b,v,i,j-1)
    else:
        return LCS(b,v,i-1,j)
    
def Alignment(b,v,w,i,j):
    # Alignment of v and w
    if ((i == 0) and (j == 0)):
        return ['','']
    if (b[i,j] == 3):
        result = Alignment(b,v,w,i-1,j-1)
        result[0] += v[i-1]
        result[1] += w[j-1]
        return result
    if (b[i,j] == 2):
        result = Alignment(b,v,w,i,j-1)
        result[0] += "_"
        result[1] += w[j-1]
        return result
    if (b[i,j] == 1):
        result = Alignment(b,v,w,i-1,j)
        result[0] += v[i-1]
        result[1] += "_"
        return result
    
def GlobalAlign(v, w, scorematrix, indel):
    # Alignment considering every character in both strings.
    s = np.zeros((len(v)+1,len(w)+1), dtype="int32")
    b = np.zeros((len(v)+1,len(w)+1), dtype="int32")
    for i in range(0,len(v)+1):
        for j in range(0,len(w)+1):
            if (j == 0):
                if (i > 0):
                    s[i,j] = s[i-1,j] + indel
                    b[i,j] = 1
                continue
            if (i == 0):
                s[i,j] = s[i,j-1] + indel
                b[i,j] = 2
                continue
            score = s[i-1,j-1] + scorematrix[v[i-1],w[j-1]]
            vskip = s[i-1,j] + indel
            wskip = s[i,j-1] + indel
            s[i,j] = max(vskip, wskip, score)
            if (s[i,j] == vskip):
                b[i,j] = 1
            elif (s[i,j] == wskip):
                b[i,j] = 2
            else:
                b[i,j] = 3
    return (s, b)

def LocalAlign(v, w, scorematrix, indel):
    s = np.zeros((len(v)+1,len(w)+1), dtype="int32")
    b = np.zeros((len(v)+1,len(w)+1), dtype="int32")
    for i in range(1,len(v)+1):
        for j in range(1,len(w)+1):
            if (j == 0):
                if (i > 0):
                    s[i,j] = max(s[i-1,j] + indel, 0)
                    b[i,j] = 1
                continue
            if (i == 0):
                s[i,j] = max(s[i,j-1] + indel, 0)
                b[i,j] = 2
                continue
            score = s[i-1,j-1] + scorematrix[v[i-1],w[j-1]]
            vskip = s[i-1,j] + indel
            wskip = s[i,j-1] + indel
            s[i,j] = max(vskip, wskip, score, 0)
            if (s[i,j] == vskip):
                b[i,j] = 1
            elif (s[i,j] == wskip):
                b[i,j] = 2
            elif (s[i,j] == score):
                b[i,j] = 3
            else:
                b[i,j] = 0
    return (s, b)

def LocalAlignment(b,v,w,i,j):
    if (b[i,j] == 0):
        return ['','']
    if (b[i,j] == 3):
        result = Alignment(b,v,w,i-1,j-1)
        result[0] += v[i-1]
        result[1] += w[j-1]
        return result
    if (b[i,j] == 2):
        result = Alignment(b,v,w,i,j-1)
        result[0] += "_"
        result[1] += w[j-1]
        return result
    if (b[i,j] == 1):
        result = Alignment(b,v,w,i-1,j)
        result[0] += v[i-1]
        result[1] += "_"
        return result