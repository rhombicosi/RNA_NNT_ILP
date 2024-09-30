import gurobipy as gp
from gurobipy import GRB
from ilp_parameters_grb import * 

# canonical base pairs
cbp_list = ['AU','UA','CG','GC','GU','UG']

def legal(RNA, i,j):
    if j - i > minD:
        if RNA[i-1] + RNA[j-1] in cbp_list:
            return 1
    else:
        return 0
    
def add_binary_vars(RNA, mip):
    listP = {}      # listP are the canonical base pairs variables
    listQ = {}      # listQ are the stacking quartets variables
    listF = {}      # listF are the first stacking quartets variables
    listL = {}      # listL are the last stacking quartets variables
    listH = {}      # listH are the hairpin loops variables
    listI = {}      # listI are the internal loops variables
    listB = {}      # listB are the bulge loop variables
    listM = {}      # listM are the 3 multi loop variables
    listX = {}      # listX are nucleotides of the hairpin loop variables
    listY = {}      # listY are nucleotides of the internal loop variables
    listZ = {}      # listZ are nucleotides of the bulge loop variables
    listW = {}      # listW are nucleotides of the bulge loop variables

    n = len(RNA)

    try:
        # create P variables
        for i in range(1, n): 
            for j in range(i + minD + 1, n + 1):
                if legal(RNA,i,j):
                        listP[f'P({i},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'P({i},{j})')                
                        
                        # listF[f'F({i},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'F({i},{j})')
        
        # create Q and F variables
        for i in range(1, n): 
            for j in range(i + minD + 1, n + 1):
                if legal(RNA,i,j):
                    if legal(RNA,i + 1, j - 1):
                            listQ[f'Q({i},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'Q({i},{j})')
                            listF[f'F({i},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'F({i},{j})')
        
        # create L variables    
        for i in range(1, n): 
            for j in range(i + minD + 1, n + 1):
                if legal(RNA,i,j):
                    if legal(RNA,i + 1, j - 1):
                        if i > 1 and j < n:
                            listL[f'L({i},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'L({i},{j})')

        # create X variables for each nucleotide
        for i in range(1, n):
            listX[f'X({i})'] = mip.addVar(vtype=GRB.BINARY, name=f'X({i})')

        # create Y variables for each nucleotide
        for i in range(1, n):
            listY[f'Y({i})'] = mip.addVar(vtype=GRB.BINARY, name=f'Y({i})')

        # create Z variables for each nucleotide
        for i in range(1, n):
            listZ[f'Z({i})'] = mip.addVar(vtype=GRB.BINARY, name=f'Z({i})')

        # create W variables for each nucleotide
        for i in range(1, n):
            listW[f'W({i})'] = mip.addVar(vtype=GRB.BINARY, name=f'W({i})')

        # create H variables
        for i in range(1,n - minD - 1):
            for j in range(i + minD + 1, n + 1):
                if RNA[i-1] + RNA[j-1] in cbp_list:
                    listH[f'H({i},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'H({i},{j})')
        
        # create I variables
        for i in range(1, n - minI - 1 - minD - 1 - minI - 1):
            for k in range(i + minI + 1, n - minI - 1 - minD - 1):
                for l in range(k + minD + 1, n - minI  - 1):
                    for j in range(l + minI + 1, n + 1):                        
                        if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                            listI[f'I({i},{k},{l},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'I({i},{k},{l},{j})')

        # create B variabels
        for i in range(1, n - 1):
            for k in range(1, n - 1):
                if k == i+1:
                    for l in range(k + minD + 1, n - 2):
                        for j in range(l + 2, n - 1):
                            if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                                listB[f'B({i},{k},{l},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'B({i},{k},{l},{j})')
                elif k == i-1:
                    for l in range(k-minD-1,4,-1):
                        for j in range(l - 2, 1, -1):
                            if RNA[j-1] + RNA[i-1] in cbp_list and RNA[l-1] + RNA[k-1] in cbp_list:
                                listB[f'B({j},{l},{k},{i})'] = mip.addVar(vtype=GRB.BINARY, name=f'B({j},{l},{k},{i})')

        # create M variables
        for i in range(1, n - 6):
            for i1 in range(i + 1, n - 5):
                for j1 in range(i1 + minD + 1, n - 4):
                    for i2 in range(j1 + 1, n - 3):
                        for j2 in range(i2 + minD + 1, n - 2):
                            for j in range(j2 + 1, n - 1):

                                # check whether all base pairs are either canonical or wobble
                                # if legal(RNA,i,j) and legal(RNA,i1,j1) and legal(RNA,i2,j2):
                                if RNA[i-1] + RNA[j-1] in cbp_list and \
                                RNA[i1-1] + RNA[j1-1] in cbp_list and \
                                RNA[i2-1] + RNA[j2-1] in cbp_list:
                                    listM[f'M({i},{i1},{j1},{i2},{j2},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'M({i},{i1},{j1},{i2},{j2},{j})')       

        # for v in mip.getVars():
        #     print(f'{v.VarName}')

    except gp.GurobiError as e:
        print(f'Error code {e.errno}: {e}')

    except AttributeError:
        print('Encountered an attribute error')  
    
    mip.update() 
    print(RNA)
    print(n)
    # print(listM)
    # print(len(listM))

    return listP, listQ, listF, listL, listH, listI, listB, listM, listX, listY, listZ, listW
