from binary_variables_grb import *
from tnn_grb import *

def num_to_range(num, inMin, inMax, outMin, outMax):
  return outMin + (float(num - inMin) / float(inMax - inMin) * (outMax
                  - outMin))

def stemParams(RNA):
    g = {}
    n = len(RNA)
    for i in range(1, n):  
        for j in range(i + minD + 1, n + 1):
            if legal(RNA,i,j):
                if legal(RNA,i + 1, j - 1):
                    g[f'Q{i}',f'Q{j}'] = G_stem(RNA,i,j)
                    
    return g

def fParams(RNA):
    g = {}
    n = len(RNA)
    for i in range(1, n):  
        for j in range(i + minD + 1, n + 1): 
            if legal(RNA,i,j):
                if legal(RNA,i+1,j-1) and i+1 < j-1:
                    g[f'F{i}',f'F{j}'] = G_F(RNA,i,j)

    return g

def lParams(RNA):
    g = {}
    n = len(RNA)
    for i in range(1, n):  
        for j in range(i + minD + 1, n + 1): 
            if legal(RNA,i,j):
                if legal(RNA,i+1,j-1) and i+1 < j-1:
                    if j < len(RNA):
                        g[f'L{i}',f'L{j}'] = G_L(RNA,i,j)

    return g 

def hairpinParams(RNA):
    g = {}
    n = len(RNA)
    for i in range(1,n - minD - 1):
        for j in range(i + minD + 1, n + 1):
            if RNA[i-1] + RNA[j-1] in cbp_list:  
                g[f'H{i}',f'H{j}'] = G_hairpin(RNA,i,j)
    
    return g

def internalParams(RNA):
    g = {}
    n = len(RNA)
    for i in range(1, n - minI - 1 - minD - 1 - minI - 1):
        for k in range(i + minI + 1, n - minI - 1 - minD - 1):
            for l in range(k + minD + 1, n - minI  - 1):
                for j in range(l + minI + 1, n + 1):
                    if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:

                       g[f'I{i}',f'I{k}',f'I{l}',f'I{j}'] = G_internal(RNA,i,k,l,j)

    return g

def bulgeParams(RNA):
    g = {}
    n = len(RNA)
    for i in range(1, n + 1):
        for k in range(1, n):
            if k == i+1:
                for l in range(k + minD+1, n - 1):
                    for j in range(l + 2, n + 1):
                        if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                            g[f'B{i}',f'B{k}',f'B{l}',f'B{j}'] = G_bulge(RNA,i,k,l,j)

            elif k == i-1:
                for l in range(k-minD-1,4,-1):
                    for j in range(l - 2, 1, -1):
                        if RNA[j-1] + RNA[i-1] in cbp_list and RNA[l-1] + RNA[k-1] in cbp_list:
                            g[f'B{j}',f'B{l}',f'B{k}',f'B{i}'] = G_bulge(RNA,j,l,k,i)
    
    return g

def multiParams(RNA):
    g = {}
    n = len(RNA)
    for i in range(1, n - 6):
        for i1 in range(i + 1, n - 5):
            for j1 in range(i1 + minD + 1, n - 4):
                for i2 in range(j1 + 1, n - 3):
                    for j2 in range(i2 + minD + 1, n - 2):
                        for j in range(j2 + 1, n + 1):
                            # if legal(RNA,i,j) and legal(RNA,i1,j1) and legal(RNA,i2,j2):
                            if RNA[i-1] + RNA[j-1] in cbp_list and \
                            RNA[i1-1] + RNA[j1-1] in cbp_list and \
                            RNA[i2-1] + RNA[j2-1] in cbp_list:   
                                g[f'M{i}',f'M{i1}',f'M{j1}',f'M{i2}',f'M{j2}',f'M{j}'] = G_multi(i,i1,j1,i2,j2,j)

    return g

def pretty(d, indent=0):
    for key, value in d.items():
        print('\t' * indent + str(key) + ' :: ' + str(value))
        if isinstance(value, dict):
            pretty(value, indent)

# pretty(stemParams(RNA))
# pretty(flParams(RNA))
# pretty(hairpinParams(RNA))
# pretty(internalParams(RNA))
# pretty(bulgeParams(RNA))

def stemTerm(RNA, listQ):
    objective = gp.LinExpr(list(stemParams(RNA).values()), list(listQ.values())) 
    return objective

def fTerm(RNA, listF):
    objective = gp.LinExpr(list(fParams(RNA).values()), list(listF.values()))
    return objective

def lTerm(RNA, listL):
    objective = gp.LinExpr(list(lParams(RNA).values()), list(listL.values()))
    return objective

def hairpinTerm(RNA, listH):
    objective = gp.LinExpr(list(hairpinParams(RNA).values()), list(listH.values()))
    return objective

def internalTerm(RNA, listI):
    objective = gp.LinExpr(list(internalParams(RNA).values()), list(listI.values()))
    return objective

def bulgeTerm(RNA, listB):
    objective = gp.LinExpr(list(bulgeParams(RNA).values()), list(listB.values()))
    return objective

def multiTerm(RNA, listM):
    objective = gp.LinExpr(list(multiParams(RNA).values()), list(listM.values()))
    return objective

# print(stemTerm(RNA))
# print(fTerm(RNA))
# print(lTerm(RNA))
# print(hairpinTerm(RNA))
# print(internalTerm(RNA))
# print(bulgeTerm(RNA))

def objectiveTerm(RNA, listQ, listH, listI, listB, listM):
    objective = stemTerm(RNA, listQ)
    objective.add(hairpinTerm(RNA, listH))
    objective.add(internalTerm(RNA, listI))
    objective.add(bulgeTerm(RNA, listB))
    objective.add(multiTerm(RNA, listM))

    return objective

# print(objectiveTerm(RNA))
