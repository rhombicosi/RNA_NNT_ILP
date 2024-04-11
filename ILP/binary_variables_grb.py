import gurobipy as gp
from gurobipy import GRB
from ilp_parameters import * 

# canonical base pairs
cbp_list = ['AU','UA','CG','GC','GU','UG']

listP = {}      # listP are the canonical base pairs variables
listQ = {}      # listQ are the stacking quartets variables
listF = {}      # listF are the first stacking quartets variables
listL = {}      # listL are the last stacking quartets variables
listH = {}      # listH are the hairpin loops variables
listI = {}      # listI are the internal loops variables
listB = {}      # listB are the bulge loop variables
listX = {}      # listX are nucleotides of the hairpin loop variables
listY = {}      # listY are nucleotides of the internal loop variables
listZ = {}      # listZ are nucleotides of the bulge loop variables


def legal(i,j):
    if j - i > minD:
        if RNA[i-1] + RNA[j-1] in cbp_list:
            return 1
    else:
        return 0

try:
    mip = gp.Model("MIP")
    
    # create P variables
    for i in range(1, len(RNA)): 
        for j in range(i + minD + 1, len(RNA) + 1):
            if legal(i,j):
                    listP[f'P({i},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'P({i},{j})')                
                    
                    # listF[f'F({i},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'F({i},{j})')
    
    # create Q and F variables
    for i in range(1, len(RNA)): 
        for j in range(i + minD + 1, len(RNA) + 1):
            if legal(i,j):
                if legal(i + 1, j - 1):
                        listQ[f'Q({i},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'Q({i},{j})')
                        listF[f'F({i},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'F({i},{j})')
    
    # create L variables    
    for i in range(1, len(RNA)): 
        for j in range(i + minD + 1, len(RNA) + 1):
            if legal(i,j):
                if legal(i + 1, j - 1):
                    if i > 1 and j < len(RNA):
                        listL[f'L({i},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'L({i},{j})')

    # create X variables for each nucleotide
    for i in range(1, len(RNA)):
        listX[f'X({i})'] = mip.addVar(vtype=GRB.BINARY, name=f'X({i})')

    # create Y variables for each nucleotide
    for i in range(1, len(RNA)):
        listY[f'Y({i})'] = mip.addVar(vtype=GRB.BINARY, name=f'Y({i})')

    # create Z variables for each nucleotide
    for i in range(1, len(RNA)):
        listZ[f'Z({i})'] = mip.addVar(vtype=GRB.BINARY, name=f'Z({i})')

    # create H variables
    for i in range(1,len(RNA) - minD - 1):
        for j in range(i + minD + 1, min(i + maxH + 1, len(RNA) + 1)):
            if RNA[i-1] + RNA[j-1] in cbp_list:
                listH[f'H({i},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'H({i},{j})')
    
    # create I variables
    for i in range(1, len(RNA) - minI - 1 - minD - 1 - minI - 1):
        for k in range(i + minI + 1, min(i + maxI + 1, len(RNA) - minI - 1 - minD - 1)):
            for l in range(k + minD + 1, len(RNA) - minI  - 1):
                for j in range(l + minI + 1, min(l + maxI + 1, len(RNA) + 1)):
                    if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                        listI[f'I({i},{k},{l},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'I({i},{k},{l},{j})')

    # create B variabels
    for i in range(1, len(RNA) - 1):
        for k in range(1, len(RNA) - 1):
            if k == i+1:
                for l in range(k + minD+1, len(RNA) - 2):
                    for j in range(l + 2, min(l + maxB, len(RNA) - 1)):
                        if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                            listB[f'B({i},{k},{l},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'B({i},{k},{l},{j})')
            elif k == i-1:
                for l in range(k-minD-1,4,-1):
                    for j in range(l - 2, max(l - maxB, 1), -1):
                        if RNA[j-1] + RNA[i-1] in cbp_list and RNA[l-1] + RNA[k-1] in cbp_list:
                            listB[f'B({j},{l},{k},{i})'] = mip.addVar(vtype=GRB.BINARY, name=f'B({j},{l},{k},{i})')
    mip.update()

    # for v in mip.getVars():
    #     print(f'{v.VarName}')

except gp.GurobiError as e:
    print(f'Error code {e.errno}: {e}')

except AttributeError:
    print('Encountered an attribute error')  

print(RNA)
print(len(RNA))




