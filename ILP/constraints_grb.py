from binary_variables_grb import *

############## BASE PAIRS :: SET OF CONSTRAINTS #############

def onePairConstraints(RNA): 
    for i in range(1,len(RNA)+1):
        inequality = gp.LinExpr(0)

        inq1 = False
        inq2 = False

        for j in range(1, i):
            if i - j > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list: 
                    inequality.add(gp.LinExpr([1.0], [mip.getVarByName(f'P({j},{i})')]))
                    inq1 = True        

        for j in range(i+1, len(RNA)+1):
            if j - i > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list:
                    inequality.add(gp.LinExpr([1.0], [mip.getVarByName(f'P({i},{j})')]))
                    inq2=True

        if inq1 or inq2:
            mip.addConstr(inequality <= 1, f'CP{i}')

    return inequality

def noCrossConstraints(RNA):  

    for h in range(1, len(RNA) - 2):
        for i in range(h+1, len(RNA) - 1):
            for j in range(i+1, len(RNA)):
                for k in range(j+1, len(RNA) + 1):
                    if j - h > minD and k - i > minD:
                        if RNA[h-1] + RNA[j-1] in cbp_list and RNA[i-1] + RNA[k-1] in cbp_list:
                            inequality = gp.LinExpr(0)
                            inequality.add(gp.LinExpr([1.0,1.0],[mip.getVarByName(f'P({h},{j})'), mip.getVarByName(f'P({i},{k})')]))
                            mip.addConstr(inequality <= 1, f'CNC{h}-{i}-{j}-{k}')

    return inequality

# ############## STEM LOOP :: SET OF CONSTRAINTS #############

def stemConstraints(RNA):    

    for i in range (1, len(RNA)): 
        for j in range(i + minD + 1,  len(RNA)+1): 
            ip1 = i + 1
            jm1 = j - 1
            if legal(i,j):
                if legal(ip1,jm1) and ip1 < jm1:                     
                    inequality = gp.LinExpr([2,-1,-1],[mip.getVarByName(f'Q({i},{j})'),mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({ip1},{jm1})')])
                    mip.addConstr(inequality <= 0, f'CS0{i}-{j}')
                    # create the inequalities to enforce the converse, although it is
                    # redundant. You can test this fact by commenting these out. Does it make
                    # difference in the Gurobi execution times?
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

def firstPairConstraints(RNA):

    for i in range (1, len(RNA)): 
        for j in range(i + minD + 1, len(RNA) + 1):            
            if legal(i,j):
                if legal(i+1,j-1):
                    if i > 1 and j < len(RNA):            
                        if legal(i-1,j+1):
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

def lastPairConstraints(RNA):

    for i in range (1, len(RNA)): 
        for j in range(i + minD + 1, len(RNA) + 1):            
            if legal(i,j):
                if legal(i+1,j-1): 
                    if legal(i+2,j-2):
                        if i > 1 and j < len(RNA):                   
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

def hairpinNTConstraints(RNA):

    for i in range(1,len(RNA)):
        inequality = gp.LinExpr(0)

        inq1 = False
        inq2 = False

        for j in range(1, i):
            if i - j > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list: 

                    inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'P({j},{i})')]))
                    inq1 = True        

        for j in range(i+1, len(RNA)):
            if j - i > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list:
                    inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'P({i},{j})')]))
                    inq2=True

        if inq1 or inq2:
            inequality.add(gp.LinExpr([1.0],[mip.getVarByName(f'X({i})')]))
            mip.addConstr(inequality >= 1, f'CX{i}')

    return inequality

def hairpinIfThenConstraints(RNA):
    
    for i in range(1,len(RNA) - minD):
        for j in range(i + minD + 1, min(i + maxH + 1, len(RNA))):
            if RNA[i-1] + RNA[j-1] in cbp_list:

                inequality = gp.LinExpr(0)
                for u in range(i+1,j):
                    inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))

                inequality.add(gp.LinExpr([1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'H({i},{j})')]))
                mip.addConstr(inequality <= j - i - 1, f'CHIFT{i}-{j}')

    return inequality

def hairpinOnlyIfConstraints(RNA):
    
    for i in range(1,len(RNA) - minD - 1):
        for j in range(i + minD + 1, min(i + maxH + 1, len(RNA) + 1)):
            if RNA[i-1] + RNA[j-1] in cbp_list:
                
                for u in range(i+1,j):
                    inequality = gp.LinExpr(0)
                    inequality.add(gp.LinExpr([2],[mip.getVarByName(f'H({i},{j})')]))
                    set1 = set({})
                    
                    for v in range(1,u):
                        if legal(v,u):
                            set1.add(f'P({v},{u})')                   

                    for v in range(u+1,len(RNA)+1):
                        if legal(u,v):
                            set1.add(f'P({u},{v})')

                    for s in set1:
                        inequality.add(gp.LinExpr([1],[mip.getVarByName(s)]))
                    
                    inequality.add(gp.LinExpr([-1], [mip.getVarByName(f'P({i},{j})')]))
                    mip.addConstr(inequality <= 1, f'CHOIF{i}-{j}-{u}')

    return inequality

def numHairpinConstraints(RNA, numH):
    
    inequality = gp.LinExpr(0)
    for i in range(1,len(RNA) - minD - 1):
        for j in range(i + minD + 1, min(i + maxH + 1, len(RNA) + 1)):
            if RNA[i-1] + RNA[j-1] in cbp_list:
                inequality.add(gp.LinExpr([1],[mip.getVarByName(f'H({i},{j})')]))

                # inequality += f' + H({i},{j})'
    mip.addConstr(inequality <= numH, f'CHN')
    # inequality += f' <= {numH} \n'

    return inequality

# ############# INTERNAL LOOP :: SET OF CONSTRAINTS ############

# unpaired nucleotides constraints
def internalNTConstraints(RNA):

    for i in range(1,len(RNA)):

        inequality = gp.LinExpr(0)

        inq1 = False
        inq2 = False

        for j in range(1, i):
            if i - j > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list: 
                    inequality.add(gp.LinExpr([1],[mip.getVarByName(f'P({j},{i})')]))
                    # inequality += f' + P({j},{i})'
                    inq1 = True        

        for j in range(i+1, len(RNA)+1):
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

# if n1 nucleotides and n2 nucleotides are unpaired between two base pairs 
# then internal loop is formed between these pairs
def internalIfThenConstraints(RNA):

    for i in range(1, len(RNA) - minI - 1 - minD - 1 - minI - 1):
        for k in range(i + minI + 1, min(i + maxI + 1, len(RNA) - minI - 1 - minD - 1)):
            for l in range(k + minD + 1, len(RNA) - minI  - 1):
                for j in range(l + minI + 1, min(l + maxI + 1, len(RNA) + 1)):
                    if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:

                        inequality = gp.LinExpr(0)                        

                        for u in range(i+1,k):
                            inequality.add(gp.LinExpr([1],[mip.getVarByName(f'Y({u})')]))

                        for u in range(l+1,j):
                            inequality.add(gp.LinExpr([1],[mip.getVarByName(f'Y({u})')]))
                            
                        inequality.add(gp.LinExpr([1,1,-1],[mip.getVarByName(f'P({k},{l})'),mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'I({i},{k},{l},{j})')]))

                        mip.addConstr(inequality <= k-i+j-l-1, f'CIIFT{i}-{k}-{l}-{j}')

    return inequality

def internalOnlyIfConstraints(RNA):
    # inequality = ''
    
    for i in range(1, len(RNA) - minI - 1 - minD - 1 - minI - 1):
        for k in range(i + minI + 1, min(i + maxI + 1, len(RNA) - minI - 1 - minD - 1)):
            for l in range(k + minD + 1, len(RNA) - minI  - 1):
                for j in range(l + minI + 1, min(l + maxI + 1, len(RNA) + 1)):
                    if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:                        

                        for u in range(i+1,k):
                            inequality = gp.LinExpr(0)

                            inequality.add(gp.LinExpr([3],[mip.getVarByName(f'I({i},{k},{l},{j})')]))
                            # inequality += f'3 I({i},{k},{l},{j}) '
                            set1 = set({})

                            for v in range(1,u):
                                if legal(v,u):
                                    set1.add(f'P({v},{u})')                   

                            for v in range(u+1,len(RNA)+1):
                                if legal(u,v):
                                    set1.add(f'P({u},{v})')

                            for s in set1:
                                inequality.add(gp.LinExpr([1],[mip.getVarByName(s)]))
                                # inequality += f'+ {s} '

                            inequality.add(gp.LinExpr([-1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({k},{l})')]))
                            mip.addConstr(inequality <= 1, f'CIOIF{i}-{k}-{l}-{j}-{u}')
                            # inequality += f'- P({i},{j}) - P({k},{l}) <= 1 \n'

                        for u in range(l+1,j):
                            inequality = gp.LinExpr(0)

                            inequality.add(gp.LinExpr([3],[mip.getVarByName(f'I({i},{k},{l},{j})')]))
                            # inequality += f'3 I({i},{k},{l},{j}) '
                            set1 = set({})

                            for v in range(1,len(RNA)+1):
                                for v in range(1,u):
                                    if legal(v,u):
                                        set1.add(f'P({v},{u})')                   

                            for v in range(u+1,len(RNA)+1):
                                if legal(u,v):
                                    set1.add(f'P({u},{v})')

                            for s in set1:
                                inequality.add(gp.LinExpr([1],[mip.getVarByName(s)]))
                                # inequality += f'+ {s} '

                            inequality.add(gp.LinExpr([-1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({k},{l})')]))
                        mip.addConstr(inequality <= 1, f'CIOIF{i}-{k}-{l}-{j}-{u}')
                            # inequality += f'- P({i},{j}) - P({k},{l}) <= 1 \n'
    
    return inequality

def numInternalConstraints(RNA, numI):
    inequality = gp.LinExpr(0)
    
    for i in range(1, len(RNA) - minI - 1 - minD - 1 - minI - 1):
        for k in range(i + minI + 1, min(i + maxI + 1, len(RNA) - minI - 1 - minD - 1)):
            for l in range(k + minD + 1, len(RNA) - minI  - 1):
                for j in range(l + minI + 1, min(l + maxI + 1, len(RNA) + 1)):
                    if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:                        
                        inequality.add(gp.LinExpr([1],[mip.getVarByName(f'I({i},{k},{l},{j})')]))

    mip.addConstr(inequality <= numI, 'CIN')

    return inequality

# ############## BULGE LOOP :: SET OF CONSTRAINTS #############

# unpaired nucleotides constraints
# TODO: repeats X and Y 
# should be joined into one function for all three variables
def bulgeNTConstraints(RNA): 

    for i in range(1,len(RNA)):

        inequality = gp.LinExpr(0)

        inq1 = False
        inq2 = False

        for j in range(1, i):
            if i - j > minD:
                if RNA[i-1] + RNA[j-1] in cbp_list: 
                    inequality.add(gp.LinExpr([1],[mip.getVarByName(f'P({j},{i})')]))
                    # inequality += f' + P({j},{i})'
                    inq1 = True        

        for j in range(i+1, len(RNA)+1):
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

# if base pairs are formed by i and j, l and k nucleotides
# and all nucleotides are unpaired between l and j
# then bulge loop is formed 
def bulgeIfThenConstraints(RNA):

    for i in range(1, len(RNA) - 1):
        for k in range(1, len(RNA) - 1):
            if k == i+1:
                for l in range(k + minD+1, len(RNA) - 2):
                    for j in range(l + 2, min(l + maxB, len(RNA) - 1)):
                        if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                            inequality = gp.LinExpr(0)                        

                            for u in range(l+1,j):
                                inequality.add(gp.LinExpr([1],[mip.getVarByName(f'Z({u})')]))

                            inequality.add(gp.LinExpr([1,1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({k},{l})'),mip.getVarByName(f'B({i},{k},{l},{j})')]))
                            mip.addConstr(inequality <= j-l, f'CBIFT{i}-{k}-{l}-{j}')

            elif k == i-1:
                for l in range(k-minD-1,4,-1):
                    for j in range(l - 2, max(l - maxB, 1), -1):
                        if RNA[j-1] + RNA[i-1] in cbp_list and RNA[l-1] + RNA[k-1] in cbp_list:
                            inequality = gp.LinExpr(0)                       

                            for u in range(j+1,l):
                                inequality.add(gp.LinExpr([1],[mip.getVarByName(f'Z({u})')]))

                            inequality.add(gp.LinExpr([1,1,-1],[mip.getVarByName(f'P({j},{i})'),mip.getVarByName(f'P({l},{k})'),mip.getVarByName(f'B({j},{l},{k},{i})')]))
                            mip.addConstr(inequality <= l-j, f'CBIFT{j}-{l}-{k}-{i}')
    
    return inequality

def bulgeOnlyIfConstraints(RNA):

    for i in range(1, len(RNA) - 1):
        for k in range(1, len(RNA) - 1):
            if k == i+1:
                for l in range(k + minD+1, len(RNA) - 2):
                    for j in range(l + 2, min(l + maxB, len(RNA) - 1)):
                        if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                            
                            for u in range(l+1,j):
                                inequality = gp.LinExpr(0)
                                inequality.add(gp.LinExpr([3],[mip.getVarByName(f'B({i},{k},{l},{j})')]))

                                set1 = set({})

                                for v in range(1,u):
                                    if legal(v,u):
                                        set1.add(f'P({v},{u})')                   

                                for v in range(u+1,len(RNA)+1):
                                    if legal(u,v):
                                        set1.add(f'P({u},{v})')

                                for s in set1:
                                    inequality.add(gp.LinExpr([1],[mip.getVarByName(s)]))

                                inequality.add(gp.LinExpr([-1,-1],[mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'P({k},{l})')]))
                                mip.addConstr(inequality <= 1, f'CBOIF{i}-{k}-{l}-{j}-{u}')

            elif k == i-1:
                for l in range(k-minD-1,4,-1):
                    for j in range(l - 2, max(l - maxB, 1), -1):
                        if RNA[j-1] + RNA[i-1] in cbp_list and RNA[l-1] + RNA[k-1] in cbp_list:                            
                            
                            for u in range(j+1,l):
                                inequality = gp.LinExpr(0)
                                inequality.add(gp.LinExpr([3],[mip.getVarByName(f'B({j},{l},{k},{i})')]))

                                set1 = set({})

                                for v in range(1,u):
                                    if legal(v,u):
                                        set1.add(f'P({v},{u})')                   

                                for v in range(u+1,len(RNA)+1):
                                    if legal(u,v):
                                        set1.add(f'P({u},{v})')

                                for s in set1:
                                    inequality.add(gp.LinExpr([1],[mip.getVarByName(s)]))

                                inequality.add(gp.LinExpr([-1,-1],[mip.getVarByName(f'P({j},{i})'),mip.getVarByName(f'P({l},{k})')]))
                                mip.addConstr(inequality <= 1, f'CBOIF{j}-{l}-{k}-{i}-{u}')
    
    return inequality

def numBulgeConstraints(RNA, numB):
    inequality = gp.LinExpr(0)

    for i in range(1, len(RNA) - 1):
        for k in range(1, len(RNA) - 1):
            if k == i+1:
                for l in range(k + minD + 1, len(RNA) - 2):
                    for j in range(l + 2, min(l + maxB, len(RNA) - 1)):
                        if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                            
                            inequality.add(gp.LinExpr([1],[mip.getVarByName(f'B({i},{k},{l},{j})')]))
                            # inequality += f' + B({i},{k},{l},{j})'

            elif k == i-1:
                for l in range(k-minD-1,4,-1):
                    for j in range(l - 2, max(l - maxB, 1), -1):
                        if RNA[j-1] + RNA[i-1] in cbp_list and RNA[l-1] + RNA[k-1] in cbp_list:
                            
                            inequality.add(gp.LinExpr([1],[mip.getVarByName(f'B({j},{l},{k},{i})')]))
                            # inequality += f' + B({j},{l},{k},{i})'

    mip.addConstr(inequality <= numB, f'CBN')
    # inequality += f' <= {numB} \n'
    
    return inequality

onePairConstraints(RNA)
noCrossConstraints(RNA)
stemConstraints(RNA)
firstPairConstraints(RNA)
lastPairConstraints(RNA)
hairpinNTConstraints(RNA)
hairpinIfThenConstraints(RNA)
hairpinOnlyIfConstraints(RNA)
numHairpinConstraints(RNA, numH)
internalNTConstraints(RNA)
internalIfThenConstraints(RNA)
internalOnlyIfConstraints(RNA)
numInternalConstraints(RNA, numI)
bulgeNTConstraints(RNA)
bulgeIfThenConstraints(RNA)
bulgeOnlyIfConstraints(RNA)
numBulgeConstraints(RNA, numB)
mip.write('1q9a_grb.lp')  

