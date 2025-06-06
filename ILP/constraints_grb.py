from binary_variables_grb import *

############## BASE PAIRS :: SET OF CONSTRAINTS #############

def onePairConstraints(RNA, mip): 
    n = len(RNA)
    for i in range(1, n + 1):
        inequality = gp.LinExpr(0)

        inq1 = False
        inq2 = False
        
        for j in range(1, i - minD):
            if i - j > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list: 
                    if j == 1 or i == n:
                        if legal(RNA,j+1,i-1):
                            inequality.add(gp.LinExpr([1.0], [mip.getVarByName(f'P({j},{i})')]))
                            inq1 = True 
                    elif legal(RNA,j+1,i-1) or legal(RNA,j-1,i+1):
                        inequality.add(gp.LinExpr([1.0], [mip.getVarByName(f'P({j},{i})')]))
                        inq1 = True 

        for j in range(i, n + 1):
            if j - i > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list:
                    if i == 1 or j == n:
                        if legal(RNA,i+1,j-1):
                            inequality.add(gp.LinExpr([1.0], [mip.getVarByName(f'P({i},{j})')]))
                            inq2=True
                    elif legal(RNA,i+1,j-1) or legal(RNA,i-1,j+1):
                        inequality.add(gp.LinExpr([1.0], [mip.getVarByName(f'P({i},{j})')]))
                        inq2=True

        if inq1 or inq2:
            mip.addConstr(inequality <= 1, f'CP{i}')

    return inequality

def noCrossConstraints(RNA, mip):
    n = len(RNA)
    for h in range(1, n - minD):
        for i in range(h+1, n - 1):
            for j in range(i+1, n):
                for k in range(j+1, n + 1):
                    if j - h > minD and k - i > minD:                        
                        if RNA[h-1] + RNA[j-1] in cbp_list and RNA[i-1] + RNA[k-1] in cbp_list:
                            inequality = gp.LinExpr(0)
                            if h == 1 or k == n:
                                if legal(RNA,h+1,j-1) and legal(RNA,i+1,k-1):                      
                                    inequality.add(gp.LinExpr([1.0,1.0],[mip.getVarByName(f'P({h},{j})'), mip.getVarByName(f'P({i},{k})')]))
                                    mip.addConstr(inequality <= 1, f'CNC{h}-{i}-{j}-{k}')
                            elif (legal(RNA,h+1,j-1) or legal(RNA,h-1,j+1)) and (legal(RNA,i+1,k-1) or legal(RNA,i-1,k+1)):
                                inequality.add(gp.LinExpr([1.0,1.0],[mip.getVarByName(f'P({h},{j})'), mip.getVarByName(f'P({i},{k})')]))
                                mip.addConstr(inequality <= 1, f'CNC{h}-{i}-{j}-{k}')

    return inequality

# ############## STEM LOOP :: SET OF CONSTRAINTS #############

def stemConstraints(RNA, mip):
    n = len(RNA)
    for i in range (1, n): 
        for j in range(i + minD + 1,  n+1): 
            ip1 = i + 1
            jm1 = j - 1
            if legal(RNA,i,j):
                if legal(RNA,ip1,jm1) and ip1 < jm1:                     
                    inequality = gp.LinExpr([2,-1,-1],[mip.getVarByName(f'Q({i},{j})'),mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({ip1},{jm1})')])
                    mip.addConstr(inequality <= 0, f'CS0{i}-{j}')
                    inequality = gp.LinExpr([1,1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({ip1},{jm1})'),mip.getVarByName(f'Q({i},{j})')])
                    mip.addConstr(inequality <= 1, f'CS1{i}-{j}')        
                # else:
                    # In this version if (i,j) is a legal pair, but (i+1,j-1) is not, the Q(i,j) is not included into the model
                    # and it is not the part of the objective function
                    # The purpose of this is to say that if (i,j) is a legal pair, but (i+1,j-1) is not, 
                    # then set Q(i,j) to 0. 
                    # This compensates for the fact that Q(i,j) may be part of the objective function.
                    # print(mip.getVarByName(f'Q({i},{j})'))
                    # inequality = gp.LinExpr([1.0],[mip.getVarByName(f'Q({i},{j})')])
                    # mip.addConstr(inequality == 0, f'CS{i}-{j}')

    return inequality

def firstPairConstraints(RNA, mip):
    n = len(RNA)
    for i in range (1, n): 
        for j in range(i + minD + 1, n + 1):            
            if legal(RNA,i,j):
                if legal(RNA,i+1,j-1):
                    if i > 1 and j < n:            
                        if legal(RNA,i-1,j+1):
                            inequality = gp.LinExpr([1,-1,-1],[mip.getVarByName(f'Q({i},{j})'),mip.getVarByName(f'Q({i-1},{j+1})'),mip.getVarByName(f'F({i},{j})')])
                            mip.addConstr(inequality <= 0, f'CF0-{i}-{j}')
                            inequality = gp.LinExpr([2,-1,1],[mip.getVarByName(f'F({i},{j})'),mip.getVarByName(f'Q({i},{j})'),mip.getVarByName(f'Q({i-1},{j+1})')])
                            mip.addConstr(inequality <= 1, f'CF1-{i}-{j}')
                    else:
                        # print(f'F: {i} :: {j}')
                        inequality = gp.LinExpr([1,-1],[mip.getVarByName(f'Q({i},{j})'),mip.getVarByName(f'F({i},{j})')])
                        mip.addConstr(inequality <= 0, f'CF0-{i}-{j}')
                        # inequality += f'Q({i},{j}) - F({i},{j}) <= 0 \n'
                        inequality = gp.LinExpr([2,-1],[mip.getVarByName(f'F({i},{j})'),mip.getVarByName(f'Q({i},{j})')])
                        mip.addConstr(inequality <= 1, f'CF1-{i}-{j}')
                        # inequality += f'2 F({i},{j}) - Q({i},{j}) <= 1 \n'

    return inequality

def lastPairConstraints(RNA, mip):
    n = len(RNA)
    for i in range (1, n): 
        for j in range(i + minD + 1, n + 1):            
            if legal(RNA,i,j):
                if legal(RNA,i+1,j-1): 
                    if legal(RNA,i+2,j-2):
                        if i > 1 and j < n:                   
                            # print(f'{i+1},{j-1} :: {mip.getVarByName(f"Q({i+1},{j-1})")}')
                            inequality = gp.LinExpr([1,-1,-1],[mip.getVarByName(f'Q({i},{j})'),mip.getVarByName(f'Q({i+1},{j-1})'),mip.getVarByName(f'L({i},{j})')])   
                            mip.addConstr(inequality <= 0, f'CL0-{i}-{j}') 
                            inequality = gp.LinExpr([2,-1,1],[mip.getVarByName(f'L({i},{j})'),mip.getVarByName(f'Q({i},{j})'),mip.getVarByName(f'Q({i+1},{j-1})')])
                            mip.addConstr(inequality <= 1, f'CL1-{i}-{j}')
                    else:
                        if i > 1 and j < len(RNA):
                            # print(f'{i},{j} :: {mip.getVarByName(f"L({i},{j})")}')
                            inequality = gp.LinExpr([1,-1],[mip.getVarByName(f'Q({i},{j})'),mip.getVarByName(f'L({i},{j})')])
                            mip.addConstr(inequality <= 0, f'CL0-{i}-{j}')
                            inequality = gp.LinExpr([2,-1],[mip.getVarByName(f'L({i},{j})'),mip.getVarByName(f'Q({i},{j})')])
                            mip.addConstr(inequality <= 1, f'CL1-{i}-{j}')

    return inequality

# ########### HAIRPIN LOOP :: SET OF CONSTRAINTS #############
def loopNTConstraints(RNA, mip):
    n = len(RNA)
    for i in range(1,n - minD):
        inequality = gp.LinExpr(0)

        inq1 = False
        inq2 = False

        for j in range(1, i):
            if i - j > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list: 
                    if j == 1:
                        if legal(RNA,j+1,i-1): # check for isolated pairs
                            inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'P({j},{i})')]))
                            inq1 = True
                    elif legal(RNA,j+1,i-1) or legal(RNA,j-1,i+1): # check for isolated pairs
                        inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'P({j},{i})')]))
                        inq1 = True     

        for j in range(i+1, n+1):
            if j - i > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list:
                    if i == 1 or j == n:
                        if legal(RNA,i+1,j-1): # check for isolated pairs
                            inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'P({i},{j})')]))
                            inq2=True
                    elif legal(RNA,i+1,j-1) or legal(RNA,i-1,j+1): # check for isolated pairs
                        inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'P({i},{j})')]))
                        inq2=True  

        if inq1 or inq2:
            inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'X({i})')]))
            mip.addConstr(inequality == 1, f'CX{i}')

    return inequality

def hairpinNTConstraints(RNA, mip): # redunduant, replaced loopNTConstraints(RNA, mip)
    n = len(RNA)
    for i in range(1,n - minD):
        inequality = gp.LinExpr(0)

        inq1 = False
        inq2 = False

        for j in range(1, i):
            if i - j > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list: 
                    if j == 1:
                        if legal(RNA,j+1,i-1): # check for isolated pairs
                            inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'P({j},{i})')]))
                            inq1 = True
                    elif legal(RNA,j+1,i-1) or legal(RNA,j-1,i+1): # check for isolated pairs
                        inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'P({j},{i})')]))
                        inq1 = True     

        for j in range(i+1, n+1):
            if j - i > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list:
                    if i == 1 or j == n:
                        if legal(RNA,i+1,j-1): # check for isolated pairs
                            inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'P({i},{j})')]))
                            inq2=True
                    elif legal(RNA,i+1,j-1) or legal(RNA,i-1,j+1): # check for isolated pairs
                        inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'P({i},{j})')]))
                        inq2=True  

        if inq1 or inq2:
            inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'X({i})')]))
            mip.addConstr(inequality == 1, f'CX{i}')

    return inequality

def hairpinZeroConstraints(RNA, mip, maxH): # redundant, replaced by consecutiveUnpairedConstraints
    n = len(RNA)
    for i in range(1, n - minD - 1):
        for j in range(i + minD + 1, n + 1):
            if RNA[i-1] + RNA[j-1] in cbp_list:                     
                if (j-i-1 > maxH):
                    inequality = gp.LinExpr([1],[mip.getVarByName(f'H({i},{j})')])
                    mip.addConstr(inequality == 0, f'HZM{i}-{j}') 

def hairpinIfThenConstraints(RNA, mip):
    n = len(RNA)
    for i in range(1,n - minD - 1):
        for j in range(i + minD + 1, n + 1):
            if RNA[i-1] + RNA[j-1] in cbp_list:
                inequality = gp.LinExpr(0)
                
                if i == 1 or j == n:
                    if legal(RNA,i+1,j-1): # check for isolated pairs
                        for u in range(i+1,j):
                            inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))

                        inequality.add(gp.LinExpr([1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'H({i},{j})')]))
                        mip.addConstr(inequality <= j - i - 1, f'CHIFT{i}-{j}')
                elif legal(RNA,i+1,j-1) or legal(RNA,i-1,j+1): # check for isolated pairs 
                    for u in range(i+1,j):
                        inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))

                    inequality.add(gp.LinExpr([1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'H({i},{j})')]))
                    mip.addConstr(inequality <= j - i - 1, f'CHIFT{i}-{j}')

    return inequality

def hairpinOnlyIfConstraints(RNA, mip): # can be ommited due to min formulation
    n = len(RNA)
    for i in range(1,n - minD - 1):
        for j in range(i + minD + 1, n + 1):
            if RNA[i-1] + RNA[j-1] in cbp_list:
                
                for u in range(i+1,j):
                    inequality = gp.LinExpr(0)
                    inequality.add(gp.LinExpr([2],[mip.getVarByName(f'H({i},{j})')]))
                    set1 = set({})
                    
                    for v in range(1,u):
                        if legal(RNA,v,u):
                            set1.add(f'P({v},{u})')                   

                    for v in range(u+1,len(RNA)+1):
                        if legal(RNA,u,v):
                            set1.add(f'P({u},{v})')

                    for s in set1:
                        inequality.add(gp.LinExpr([1],[mip.getVarByName(s)]))
                    
                    inequality.add(gp.LinExpr([-1], [mip.getVarByName(f'P({i},{j})')]))
                    mip.addConstr(inequality <= 1, f'CHOIF{i}-{j}-{u}')

    return inequality

def numHairpinConstraints(RNA, numH, mip):
    n = len(RNA)
    inequality = gp.LinExpr(0)
    for i in range(1,n - minD - 1):
        for j in range(i + minD + 1, n + 1):
            if RNA[i-1] + RNA[j-1] in cbp_list:
                inequality.add(gp.LinExpr([1],[mip.getVarByName(f'H({i},{j})')]))

    mip.addConstr(inequality <= numH, f'CHN')
    
    return inequality

# ############# INTERNAL LOOP :: SET OF CONSTRAINTS ############

# unpaired nucleotides constraints
def internalNTConstraints(RNA, mip): # redunduant, replaced loopNTConstraints(RNA, mip)
    n = len(RNA)
    for i in range(1,n):

        inequality = gp.LinExpr(0)

        inq1 = False
        inq2 = False

        for j in range(1, i):
            if i - j > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list: 
                    inequality.add(gp.LinExpr([1],[mip.getVarByName(f'P({j},{i})')]))
                    # inequality += f' + P({j},{i})'
                    inq1 = True        

        for j in range(i+1, n+1):
            if j - i > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list:
                    inequality.add(gp.LinExpr([1],[mip.getVarByName(f'P({i},{j})')]))
                    # inequality += f' + P({i},{j})'
                    inq2=True

        if inq1 or inq2:
            inequality.add(gp.LinExpr([1],[mip.getVarByName(f'Y({i})')]))
            mip.addConstr(inequality >= 1, f'CY{i}')
            # inequality += f' + Y({i}) >= 1\n'

    return inequality

def internalZeroConstraints(RNA, mip, maxI): # redundant, replaced by consecutiveUnpairedConstraints
    n = len(RNA)
    for i in range(1, n - minI - 1 - minD - 1 - minI - 1):
        for k in range(i + minI + 1, n - minI - 1 - minD - 1):
            for l in range(k + minD + 1, n - minI  - 1):
                for j in range(l + minI + 1, n + 1):
                    if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                        
                        if (k-i-1 > maxI) or (j-l-1 > maxI):
                                inequality = gp.LinExpr([1],[mip.getVarByName(f'I({i},{k},{l},{j})')])
                                mip.addConstr(inequality == 0, f'IZM{i}-{k}-{l}-{j}') 

# if n1 nucleotides and n2 nucleotides are unpaired between two base pairs 
# then internal loop is formed between these pairs
def internalIfThenConstraints(RNA, mip):
    n = len(RNA)
    for i in range(1, n - minI - 1 - minD - 1 - minI - 1):
        for k in range(i + minI + 1, n - minI - 1 - minD - 1):
            for l in range(k + minD + 1, n - minI  - 1):
                for j in range(l + minI + 1, n + 1):
                    if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                        inequality = gp.LinExpr(0)                        

                        if i == 1 or j == n:
                            pass
                        elif legal(RNA,k+1,l-1) and legal(RNA,i-1,j+1): # check for isolated pairs 

                            for u in range(i+1,k):
                                inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))

                            for u in range(l+1,j):
                                inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))
                                
                            inequality.add(gp.LinExpr([1,1,-1],[mip.getVarByName(f'P({k},{l})'),mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'I({i},{k},{l},{j})')]))

                            mip.addConstr(inequality <= round(k-i+j-l-1,2), f'CIIFT{i}-{k}-{l}-{j}')                      

                            # for u in range(i+1,k):
                            #     inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))

                            # for u in range(l+1,j):
                            #     inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))
                                
                            # inequality.add(gp.LinExpr([1,1,-1],[mip.getVarByName(f'P({k},{l})'),mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'I({i},{k},{l},{j})')]))

                            # mip.addConstr(inequality <= round(k-i+j-l-1,2), f'CIIFT{i}-{k}-{l}-{j}')

    return inequality

def internalOnlyIfConstraints(RNA, mip):
    # inequality = ''
    n = len(RNA)
    for i in range(1, n - minI - 1 - minD - 1 - minI - 1):
        for k in range(i + minI + 1, n - minI - 1 - minD - 1):
            for l in range(k + minD + 1, n - minI  - 1):
                for j in range(l + minI + 1, n + 1):
                    if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:

                        if i == 1 or j == n:
                            pass
                        elif legal(RNA,k+1,l-1) and legal(RNA,i-1,j+1): # check for isolated pairs
                            for u in range(i+1,k):
                                inequality = gp.LinExpr(0)
                                inequality.add(gp.LinExpr([3],[mip.getVarByName(f'I({i},{k},{l},{j})')]))                                
                                set1 = set({})

                                for v in range(1,u):
                                    if legal(RNA,v,u):
                                        if v == 1:
                                            if legal(RNA,v+1,u-1): # check for isolated pairs   
                                                set1.add(f'P({v},{u})')
                                        elif legal(RNA,v+1,u-1) or legal(RNA,v-1,u+1): # check for isolated pairs
                                            set1.add(f'P({v},{u})')                   

                                for v in range(u+1,n+1):
                                    if legal(RNA,u,v):
                                        if v == n:
                                            if legal(RNA,u+1,v-1): # check for isolated pairs
                                                set1.add(f'P({u},{v})')
                                        elif legal(RNA,u+1,v-1) or legal(RNA,u-1,v+1): # check for isolated pairs
                                            set1.add(f'P({u},{v})')

                                for s in set1:
                                    inequality.add(gp.LinExpr([1],[mip.getVarByName(s)]))

                                inequality.add(gp.LinExpr([-1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({k},{l})')]))
                                mip.addConstr(inequality <= 1, f'CIOIF{i}-{k}-{l}-{j}-{u}')

                            for u in range(l+1,j):
                                inequality = gp.LinExpr(0)
                                inequality.add(gp.LinExpr([3],[mip.getVarByName(f'I({i},{k},{l},{j})')]))                                
                                set1 = set({})

                                if v == 1:
                                    if legal(RNA,v+1,u-1): # check for isolated pairs   
                                        set1.add(f'P({v},{u})')
                                elif legal(RNA,v+1,u-1) or legal(RNA,v-1,u+1): # check for isolated pairs
                                    set1.add(f'P({v},{u})') 
                                
                                for v in range(u+1,n+1):
                                    if legal(RNA,u,v):
                                        if v == n:
                                            if legal(RNA,u+1,v-1): # check for isolated pairs
                                                set1.add(f'P({u},{v})')
                                        elif legal(RNA,u+1,v-1) or legal(RNA,u-1,v+1): # check for isolated pairs
                                            set1.add(f'P({u},{v})')

                                for s in set1:
                                    inequality.add(gp.LinExpr([1],[mip.getVarByName(s)]))

                                inequality.add(gp.LinExpr([-1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({k},{l})')]))
                                mip.addConstr(inequality <= 1, f'CIOIF{i}-{k}-{l}-{j}-{u}')
    
    return inequality

def numInternalConstraints(RNA, numI, mip):
    inequality = gp.LinExpr(0)
    n = len(RNA)
    for i in range(1, n - minI - 1 - minD - 1 - minI - 1):
        for k in range(i + minI + 1, n - minI - 1 - minD - 1):
            for l in range(k + minD + 1, n - minI  - 1):
                for j in range(l + minI + 1, n + 1):
                    if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:                        
                        inequality.add(gp.LinExpr([1],[mip.getVarByName(f'I({i},{k},{l},{j})')]))

    mip.addConstr(inequality >= numI, 'CIN')

    return inequality

# ############## BULGE LOOP :: SET OF CONSTRAINTS #############

# unpaired nucleotides constraints
# TODO: repeats X and Y 
# should be joined into one function for all three variables
def bulgeNTConstraints(RNA, mip): # redunduant, replaced loopNTConstraints(RNA, mip) 
    n = len(RNA)
    for i in range(1,n):

        inequality = gp.LinExpr(0)

        inq1 = False
        inq2 = False

        for j in range(1, i):
            if i - j > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list: 
                    inequality.add(gp.LinExpr([1],[mip.getVarByName(f'P({j},{i})')]))
                    # inequality += f' + P({j},{i})'
                    inq1 = True        

        for j in range(i+1, n+1):
            if j - i > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list:
                    inequality.add(gp.LinExpr([1],[mip.getVarByName(f'P({i},{j})')]))
                    # inequality += f' + P({i},{j})'
                    inq2=True

        if inq1 or inq2:
            inequality.add(gp.LinExpr([1],[mip.getVarByName(f'Z({i})')]))
            mip.addConstr(inequality >= 1, f'CZ{i}')
            # inequality += f' + Z({i}) >= 1\n'

    return inequality

# zero all variable that correspond to loops bigger then noB and maxB
def bulgeZeroConstraints(RNA, mip, maxB): # redundant, replaced by consecutiveUnpairedConstraints
    n = len(RNA)
    for i in range(1, n+1):
        for k in range(1, n):
            if k == i+1:
                for l in range(k + minD + 1, n - 1):
                    for j in range(l + 2, n + 1):
                        if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                            if (j-l-1 > maxB):
                                inequality = gp.LinExpr([1],[mip.getVarByName(f'B({i},{k},{l},{j})')])                                
                                mip.addConstr(inequality == 0, f'CZB{i}-{k}-{l}-{j}')
            elif k == i-1:
                for l in range(k-minD-1,4,-1):
                    for j in range(l - 2, 1, -1):
                        if RNA[j-1] + RNA[i-1] in cbp_list and RNA[l-1] + RNA[k-1] in cbp_list:
                            if (l-j-1 > maxB):
                                inequality = gp.LinExpr([1],[mip.getVarByName(f'B({j},{l},{k},{i})')])
                                mip.addConstr(inequality == 0, f'CZB{j}-{l}-{k}-{i}')

# if base pairs are formed by i and j, l and k nucleotides
# and all nucleotides are unpaired between l and j
# then bulge loop is formed 
def bulgeIfThenConstraints(RNA, mip):
    n = len(RNA)
    for i in range(1, n - minD + 1):
        for k in range(1, n):
            if k == i+1:
                for l in range(k + minD+1, n-1):
                    for j in range(l + 2, n+1):
                        if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:                            
                            inequality = gp.LinExpr(0)                                                

                            if i == 1 or j == n:
                                pass
                            elif legal(RNA,k+1,l-1) and legal(RNA,i-1,j+1): # check for isolated pairs
                                for u in range(l+1,j):
                                    inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))

                                inequality.add(gp.LinExpr([1,1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({k},{l})'),mip.getVarByName(f'B({i},{k},{l},{j})')]))
                                mip.addConstr(inequality <= j-l, f'CBIFT{i}-{k}-{l}-{j}')
                            
                            #     for u in range(l+1,j):
                            #         inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))

                            #     inequality.add(gp.LinExpr([1,1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({k},{l})'),mip.getVarByName(f'B({i},{k},{l},{j})')]))
                            #     mip.addConstr(inequality <= j-l, f'CBIFT{i}-{k}-{l}-{j}')

            elif k == i-1:
                for l in range(k-minD-1,4,-1):
                    for j in range(l - 2, 1, -1):
                        if RNA[j-1] + RNA[i-1] in cbp_list and RNA[l-1] + RNA[k-1] in cbp_list:
                            inequality = gp.LinExpr(0)

                            if j == 1 or i == n:
                                pass
                            elif legal(RNA,l+1,k-1) and legal(RNA,j-1,i+1): # check for isolated pairs
                                for u in range(j+1,l):
                                    inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))

                                inequality.add(gp.LinExpr([1,1,-1],[mip.getVarByName(f'P({j},{i})'),mip.getVarByName(f'P({l},{k})'),mip.getVarByName(f'B({j},{l},{k},{i})')]))
                                mip.addConstr(inequality <= l-j, f'CBIFT{j}-{l}-{k}-{i}')
    
    return inequality

def bulgeOnlyIfConstraints(RNA, mip):
    n = len(RNA)
    for i in range(1, n - minD + 1):
        for k in range(1, n):
            if k == i+1:
                for l in range(k + minD+1, n - 1):
                    for j in range(l + 2, n + 1):
                        if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                            
                            for u in range(l+1,j):
                                inequality = gp.LinExpr(0)
                                inequality.add(gp.LinExpr([3],[mip.getVarByName(f'B({i},{k},{l},{j})')]))

                                set1 = set({})

                                for v in range(1,u):
                                    if legal(RNA,v,u):
                                        set1.add(f'P({v},{u})')                   

                                for v in range(u+1,n+1):
                                    if legal(RNA,u,v):
                                        set1.add(f'P({u},{v})')

                                for s in set1:
                                    inequality.add(gp.LinExpr([1],[mip.getVarByName(s)]))

                                inequality.add(gp.LinExpr([-1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({k},{l})')]))
                                mip.addConstr(inequality <= 1, f'CBOIF{i}-{k}-{l}-{j}-{u}')

            elif k == i-1:
                for l in range(k-minD-1,4,-1):
                    for j in range(l - 2, 1, -1):
                        if RNA[j-1] + RNA[i-1] in cbp_list and RNA[l-1] + RNA[k-1] in cbp_list:                            
                            
                            for u in range(j+1,l):
                                inequality = gp.LinExpr(0)
                                inequality.add(gp.LinExpr([3],[mip.getVarByName(f'B({j},{l},{k},{i})')]))

                                set1 = set({})

                                for v in range(1,u):
                                    if legal(RNA,v,u):
                                        set1.add(f'P({v},{u})')                   

                                for v in range(u+1,n+1):
                                    if legal(RNA,u,v):
                                        set1.add(f'P({u},{v})')

                                for s in set1:
                                    inequality.add(gp.LinExpr([1],[mip.getVarByName(s)]))

                                inequality.add(gp.LinExpr([-1,-1],[mip.getVarByName(f'P({j},{i})'),mip.getVarByName(f'P({l},{k})')]))
                                mip.addConstr(inequality <= 1, f'CBOIF{j}-{l}-{k}-{i}-{u}')
    
    return inequality

def numBulgeConstraints(RNA, numB, mip):
    n = len(RNA)
    inequality = gp.LinExpr(0)
    for i in range(1, n +1):
        for k in range(1, n):
            if k == i+1:
                for l in range(k + minD + 1, n - 1):
                    for j in range(l + 2, n +1):
                        if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:                            
                            inequality.add(gp.LinExpr([1],[mip.getVarByName(f'B({i},{k},{l},{j})')]))

            elif k == i-1:
                for l in range(k-minD-1,4,-1):
                    for j in range(l - 2, 1, -1):
                        if RNA[j-1] + RNA[i-1] in cbp_list and RNA[l-1] + RNA[k-1] in cbp_list:                            
                            inequality.add(gp.LinExpr([1],[mip.getVarByName(f'B({j},{l},{k},{i})')]))

    mip.addConstr(inequality <= numB, f'CBN')
    # inequality += f' <= {numB} \n'
    
    return inequality

# ############# 3 MULTI LOOP :: SET OF CONSTRAINTS ############

# unpaired nucleotides constraints
def multiNTConstraints(RNA, mip): # redunduant, replaced loopNTConstraints(RNA, mip)
    n = len(RNA)

    for i in range(1,n):

        inequality = gp.LinExpr(0)

        inq1 = False
        inq2 = False

        for j in range(1, i):
            if i - j > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list: 
                    inequality.add(gp.LinExpr([1],[mip.getVarByName(f'P({j},{i})')]))
                    inq1 = True        

        for j in range(i+1, n+1):
            if j - i > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list:
                    inequality.add(gp.LinExpr([1],[mip.getVarByName(f'P({i},{j})')]))
                    inq2=True

        if inq1 or inq2:
            inequality.add(gp.LinExpr([1],[mip.getVarByName(f'W({i})')]))
            mip.addConstr(inequality >= 1, f'CW{i}')

    return inequality

# zero all variable that correspond to loops bigger then noM and maxM
def multiZeroConstraints(RNA, mip, maxM): # redundant, replaced by consecutiveUnpairedConstraints
    n = len(RNA)

    for i in range(1, n - 6):
        for i1 in range(i + 1, n - 5):
            for j1 in range(i1 + minD + 1, n - 4):
                for i2 in range(j1 + 1, n - 3):
                    for j2 in range(i2 + minD + 1, n - 2):
                        for j in range(j2 + 1, n - 1):
                            # if legal(RNA,i,j) and legal(RNA,i1,j1) and legal(RNA,i2,j2):   
                            if RNA[i-1] + RNA[j-1] in cbp_list and \
                            RNA[i1-1] + RNA[j1-1] in cbp_list and \
                            RNA[i2-1] + RNA[j2-1] in cbp_list:                

                                if (i1-i-1 > maxM) or (i2-j1-1 > maxM) or (j2-j-1 > maxM):
                                    inequality = gp.LinExpr([1],[mip.getVarByName(f'M({i},{i1},{j1},{i2},{j2},{j})')])
                                    mip.addConstr(inequality == 0, f'CZM{i}-{i1}-{j1}-{i2}-{j2}-{j}')

    return inequality

# if n1 nucleotides and n2 nucleotides are unpaired between two base pairs 
# then multi loop is formed between these pairs
def multiIfThenConstraints(RNA, mip):
    n = len(RNA)

    for i in range(1, n - 6):
        for i1 in range(i + 1, n - 5):
            for j1 in range(i1 + minD + 1, n - 4):
                for i2 in range(j1 + 1, n - 3):
                    for j2 in range(i2 + minD + 1, n - 2):
                        for j in range(j2 + 1, n + 1):
                            if RNA[i-1] + RNA[j-1] in cbp_list and \
                            RNA[i1-1] + RNA[j1-1] in cbp_list and \
                            RNA[i2-1] + RNA[j2-1] in cbp_list:

                                inequality = gp.LinExpr(0)

                                if i == 1 or j == n:
                                    pass
                                elif legal(RNA,i-1,j+1) and legal(RNA,i1-1,j1+1) and legal(RNA,i2-1,j2+1): # check for isolated pairs
                                    for u in range(i+1,i1):
                                        inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))

                                    for u in range(j1+1,i2):
                                        inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))

                                    for u in range(j2+1,j):
                                        inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))
                                        
                                    inequality.add(gp.LinExpr([1,1,1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({i1},{j1})'),mip.getVarByName(f'P({i2},{j2})'),mip.getVarByName(f'M({i},{i1},{j1},{i2},{j2},{j})')]))

                                    mip.addConstr(inequality <= i1-i+i2-j1+j-j2-1, f'CMIFT{i}-{i1}-{j1}-{i2}-{j2}-{j}')

    return inequality

def multiOnlyIfConstraints(RNA, mip):
    n = len(RNA)

    for i in range(1, n - 6):
        for i1 in range(i + 1, n - 5):
            for j1 in range(i1 + minD + 1, n - 4):
                for i2 in range(j1 + 1, n - 3):
                    for j2 in range(i2 + minD + 1, n - 2):
                        for j in range(j2 + 1, n + 1):
                            if RNA[i-1] + RNA[j-1] in cbp_list and \
                            RNA[i1-1] + RNA[j1-1] in cbp_list and \
                            RNA[i2-1] + RNA[j2-1] in cbp_list: 


                                # if i == 1 or j == n:
                                #     if legal(RNA,i+1,j-1): # check for isolated pairs
                                #         listP[f'P({i},{j})'] = mip.addVar(vtype=GRB.BINARY, name=f'P({i},{j})')
                                # elif legal(RNA,i+1,j-1) or legal(RNA,i-1,j+1): # check for isolated pairs                      

                                for u in range(i+1,i1):
                                    inequality = gp.LinExpr(0)

                                    inequality.add(gp.LinExpr([4],[mip.getVarByName(f'M({i},{i1},{j1},{i2},{j2},{j})')]))
                                
                                    set1 = set({})

                                    for v in range(1,u):
                                        if legal(RNA,v,u):
                                            set1.add(f'P({v},{u})')                   

                                    for v in range(u+1,n+1):
                                        if legal(RNA,u,v):
                                            set1.add(f'P({u},{v})')

                                    for s in set1:
                                        inequality.add(gp.LinExpr([1],[mip.getVarByName(s)]))

                                    inequality.add(gp.LinExpr([-1,-1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({i1},{j1})'),mip.getVarByName(f'P({i2},{j2})')]))
                                    mip.addConstr(inequality <= 1, f'CMOIF{i}-{i1}-{j1}-{i2}-{j2}-{j}-{u}')

                                for u in range(j1+1,i2):
                                    inequality = gp.LinExpr(0)

                                    inequality.add(gp.LinExpr([4],[mip.getVarByName(f'M({i},{i1},{j1},{i2},{j2},{j})')]))
                                    set1 = set({})

                                    for v in range(1,n+1):
                                        for v in range(1,u):
                                            if legal(RNA,v,u):
                                                set1.add(f'P({v},{u})')                   

                                    for v in range(u+1,n+1):
                                        if legal(RNA,u,v):
                                            set1.add(f'P({u},{v})')

                                    for s in set1:
                                        inequality.add(gp.LinExpr([1],[mip.getVarByName(s)]))

                                    inequality.add(gp.LinExpr([-1,-1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({i1},{j1})'),mip.getVarByName(f'P({i2},{j2})')]))
                                    mip.addConstr(inequality <= 1, f'CMOIF{i}-{i1}-{j1}-{i2}-{j2}-{j}-{u}')

                                for u in range(j2+1,j):
                                    inequality = gp.LinExpr(0)

                                    inequality.add(gp.LinExpr([4],[mip.getVarByName(f'M({i},{i1},{j1},{i2},{j2},{j})')]))
                                    set1 = set({})

                                    for v in range(1,n+1):
                                        for v in range(1,u):
                                            if legal(RNA,v,u):
                                                set1.add(f'P({v},{u})')                   

                                    for v in range(u+1,n+1):
                                        if legal(RNA,u,v):
                                            set1.add(f'P({u},{v})')

                                    for s in set1:
                                        inequality.add(gp.LinExpr([1],[mip.getVarByName(s)]))

                                    inequality.add(gp.LinExpr([-1,-1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({i1},{j1})'),mip.getVarByName(f'P({i2},{j2})')]))
                                    mip.addConstr(inequality <= 1, f'CMOIF{i}-{i1}-{j1}-{i2}-{j2}-{j}-{u}')
    
    return inequality

def numMultiConstraints(RNA, numM, mip):
    inequality = gp.LinExpr(0)
    
    n = len(RNA)

    for i in range(1, n - 6):
        for i1 in range(i + 1, n - 5):
            for j1 in range(i1 + minD + 1, n - 4):
                for i2 in range(j1 + 1, n - 3):
                    for j2 in range(i2 + minD + 1, n - 2):
                        for j in range(j2 + 1, n - 1):
                            if RNA[i-1] + RNA[j-1] in cbp_list and \
                            RNA[i1-1] + RNA[j1-1] in cbp_list and \
                            RNA[i2-1] + RNA[j2-1] in cbp_list:                         
                                inequality.add(gp.LinExpr([1],[mip.getVarByName(f'M({i},{i1},{j1},{i2},{j2},{j})')]))

    mip.addConstr(inequality <= numM, 'CMN')

    return inequality

def numLowerMultiConstraints(RNA, numM, mip):
    inequality = gp.LinExpr(0)
    
    n = len(RNA)

    for i in range(1, n - 6):
        for i1 in range(i + 1, n - 5):
            for j1 in range(i1 + minD + 1, n - 4):
                for i2 in range(j1 + 1, n - 3):
                    for j2 in range(i2 + minD + 1, n - 2):
                        for j in range(j2 + 1, n - 1):
                            if RNA[i-1] + RNA[j-1] in cbp_list and \
                            RNA[i1-1] + RNA[j1-1] in cbp_list and \
                            RNA[i2-1] + RNA[j2-1] in cbp_list:                         
                                inequality.add(gp.LinExpr([1],[mip.getVarByName(f'M({i},{i1},{j1},{i2},{j2},{j})')]))

    mip.addConstr(inequality >= numM, 'CMGN')

    return inequality

# prohibit unpaired regions longer than L
# that is H, I, B, M loops less than L
def consecutiveUnpairedConstraints(RNA, L, mip):

    n = len(RNA)
    step = minD+2

    for i in range(1,n+1,step):
        inequality = gp.LinExpr(0)
        for k in range(i,min(i+L, n+1)):
            for j in range(1, i):
                if k - j > minD:
                    if RNA[k-1] + RNA[j-1] in cbp_list:
                        if j == 1 or k == n:
                            if legal(RNA,j+1,k-1): # check for isolated pairs 
                                inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'P({j},{k})')]))
                        elif legal(RNA,j+1,k-1) or legal(RNA,j-1,k+1): # check for isolated pairs
                            inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'P({j},{k})')]))


            for j in range(i+1, n+1):
                if j - k > minD:
                    if RNA[j-1] + RNA[k-1] in cbp_list:
                        if k == 1 or j == n:
                            if legal(RNA,k+1,j-1): # check for isolated pairs 
                                inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'P({k},{j})')]))
                        elif legal(RNA,k+1,j-1) or legal(RNA,k-1,j+1): # check for isolated pairs
                            inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'P({k},{j})')]))

        mip.addConstr(inequality >= 1, f'CUC{i}')

def isolatedPairsConstraints(RNA, mip):
    n = len(RNA)
    for i in range(1, n - minD): 
        for j in range(i + minD + 1, n + 1):
            if legal(RNA,i,j):

                inequality = gp.LinExpr(0)

                if i == 1 or j == n:
                    if legal(RNA,i+1,j-1): # check for isolated pairs
                        inequality.add(gp.LinExpr([1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({i+1},{j-1})')]))
                        mip.addConstr(inequality <= 0, f'CIP{i}-{j}')
                elif legal(RNA,i+1,j-1) and legal(RNA,i-1,j+1): # check for isolated pairs
                        inequality.add(gp.LinExpr([1,-1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({i+1},{j-1})'),mip.getVarByName(f'P({i-1},{j+1})')]))
                        mip.addConstr(inequality <= 0, f'CIP{i}-{j}')      


# : :: ::: :::: test :::: ::: :: :

# onePairConstraints(RNA)
# noCrossConstraints(RNA)
# stemConstraints(RNA)
# firstPairConstraints(RNA)
# lastPairConstraints(RNA)
# hairpinNTConstraints(RNA)
# hairpinIfThenConstraints(RNA)
# hairpinOnlyIfConstraints(RNA)
# numHairpinConstraints(RNA, numH)
# internalNTConstraints(RNA)
# internalIfThenConstraints(RNA)
# internalOnlyIfConstraints(RNA)
# numInternalConstraints(RNA, numI)
# bulgeNTConstraints(RNA)
# bulgeIfThenConstraints(RNA)
# bulgeOnlyIfConstraints(RNA)
# numBulgeConstraints(RNA, numB)


# mip.write('1q9a_grb.lp')  

