from binary_variables_grb import *
from tnn_grb import *

def num_to_range(num, inMin, inMax, outMin, outMax):
  return outMin + (float(num - inMin) / float(inMax - inMin) * (outMax
                  - outMin))

def stemParams(RNA):
    g = {}

    for i in range(1, len(RNA)):  
        for j in range(i + minD + 1, len(RNA) + 1):
            if legal(i,j):
                if legal(i + 1, j - 1):
                    g[f'Q{i}',f'Q{j}'] = G_stem(i,j)
                    
    return g

def fParams(RNA):
    g = {}

    for i in range(1, len(RNA)):  
        for j in range(i + minD + 1, len(RNA) + 1): 
            if legal(i,j):
                if legal(i+1,j-1) and i+1 < j-1:
                    g[f'F{i}',f'F{j}'] = G_F(i,j)

    return g

def lParams(RNA):
    g = {}

    for i in range(1, len(RNA)):  
        for j in range(i + minD + 1, len(RNA) + 1): 
            if legal(i,j):
                if legal(i+1,j-1) and i+1 < j-1:
                    if j < len(RNA):
                        g[f'L{i}',f'L{j}'] = G_L(i,j)

    return g 

def hairpinParams(RNA):
    g = {}

    for i in range(1,len(RNA) - minD - 1):
        for j in range(i + minD + 1, min(i + maxH + 1, len(RNA) + 1)):
            if RNA[i-1] + RNA[j-1] in cbp_list:  
                g[f'H{i}',f'H{j}'] = G_hairpin(i,j)
    
    return g

def internalParams(RNA):
    g = {}

    for i in range(1, len(RNA) - minI - 1 - minD - 1 - minI - 1):
        for k in range(i + minI + 1, min(i + maxI + 1, len(RNA) - minI - 1 - minD - 1)):
            for l in range(k + minD + 1, len(RNA) - minI  - 1):
                for j in range(l + minI + 1, min(l + maxI + 1, len(RNA) + 1)):
                    if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:

                       g[f'I{i}',f'I{k}',f'I{l}',f'I{j}'] = G_internal(i,k,l,j)

    return g

def bulgeParams(RNA):
    g = {}

    for i in range(1, len(RNA) - 1):
        for k in range(1, len(RNA) - 1):
            if k == i+1:
                for l in range(k + minD+1, len(RNA) - 2):
                    for j in range(l + 2, min(l + maxB, len(RNA) - 1)):
                        if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                            g[f'B{i}',f'B{k}',f'B{l}',f'B{j}'] = G_bulge(i,k,l,j)

            elif k == i-1:
                for l in range(k-minD-1,4,-1):
                    for j in range(l - 2, max(l - maxB, 1), -1):
                        if RNA[j-1] + RNA[i-1] in cbp_list and RNA[l-1] + RNA[k-1] in cbp_list:
                            g[f'B{j}',f'B{l}',f'B{k}',f'B{i}'] = G_bulge(j,l,k,i)
    
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

def stemTerm(RNA):
    objective = gp.LinExpr(list(stemParams(RNA).values()), list(listQ.values())) 
    return objective

def fTerm(RNA):
    objective = gp.LinExpr(list(fParams(RNA).values()), list(listF.values()))
    return objective

def lTerm(RNA):
    objective = gp.LinExpr(list(lParams(RNA).values()), list(listL.values()))
    return objective

def hairpinTerm(RNA):
    objective = gp.LinExpr(list(hairpinParams(RNA).values()), list(listH.values()))
    return objective

def internalTerm(RNA):
    objective = gp.LinExpr(list(internalParams(RNA).values()), list(listI.values()))
    return objective

def bulgeTerm(RNA):
    objective = gp.LinExpr(list(bulgeParams(RNA).values()), list(listB.values()))    
    return objective

# print(stemTerm(RNA))
# print(fTerm(RNA))
# print(lTerm(RNA))
# print(hairpinTerm(RNA))
# print(internalTerm(RNA))
# print(bulgeTerm(RNA))

def objectiveTerm(RNA):
    objective = stemTerm(RNA)
    objective.add(hairpinTerm(RNA))
    objective.add(internalTerm(RNA))
    objective.add(bulgeTerm(RNA))
    # objective = stemTerm(RNA) + hairpinTerm(RNA) + internalTerm(RNA) + bulgeTerm(RNA) #+ fTerm(RNA) + lTerm(RNA)

    return objective

# print(objectiveTerm(RNA))

try:
    mip.setObjective(objectiveTerm(RNA), GRB.MINIMIZE)

    mip.update()

    mip.optimize()

    print(f'Obj: {mip.ObjVal:g}')

except gp.GurobiError as e:
    print(f'Error code {e.errno}: {e}')

except AttributeError:
    print('Encountered an attribute error')