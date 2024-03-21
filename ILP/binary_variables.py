from ilp_parameters import * 
# import os

# read sequence
# chain_dir = 'data'
# chain_f = '2ZY6_A'
# chain_path = os.path.join(chain_dir, chain_f)

# with open( chain_path, 'r') as file:
#     RNA = file.read().rstrip()    
#     print (RNA)

# canonical base pairs
cbp_list = ['AU','UA','CG','GC','GU','UG']

listP = ""      # listP are the canonical base pairs variables
listQ = ""      # listQ are the stacking quartets variables
listF = ""      # listF are the first stacking quartets variables
listL = ""      # listL are the last stacking quartets variables
listH = ""      # listH are the hairpin loops variables
listI = ""      # listI are the internal loops variables
listB = ""      # listB are the bulge loop variables
listX = ""      # listX are nucleotides of the hairpin loop variables
listY = ""      # listY are nucleotides of the internal loop variables
listZ = ""      # listZ are nucleotides of the bulge loop variables

def legal(i,j):
    if j - i > minD:
        if RNA[i-1] + RNA[j-1] in cbp_list:
            return 1
    else:
        return 0

# create P and F variables
for i in range(1, len(RNA)): 
    for j in range(i + minD + 1, len(RNA) + 1):
        if RNA[i-1] + RNA[j-1] in cbp_list:
            listP = listP + f'P({i},{j})\n'
            listF = listF + f'F({i},{j})\n'

# create Q and L variables
for i in range(1, len(RNA)): 
    for j in range(i + minD + 1, len(RNA) + 1):
        if RNA[i-1] + RNA[j-1] and RNA[i] + RNA[j-2] in cbp_list:
            listQ = listQ + f'Q({i},{j})\n'
            listL = listL + f'L({i},{j})\n'

# create X variables for each nucleotide
for i in range(1, len(RNA)):
    listX += f'X({i})\n'

# create Y variables for each nucleotide
for i in range(1, len(RNA)):
    listY += f'Y({i})\n'

# create Z variables for each nucleotide
for i in range(1, len(RNA)):
    listZ += f'Z({i})\n'

# create H variables
for i in range(1,len(RNA) - minD - 1):
    for j in range(i + minD + 1, min(i + maxH + 1, len(RNA) + 1)):
        if RNA[i-1] + RNA[j-1] in cbp_list:
            listH += f'H({i},{j})\n'

# create I variables
for i in range(1, len(RNA) - minI - 1 - minD - 1 - minI - 1):
    for k in range(i + minI + 1, min(i + maxI + 1, len(RNA) - minI - 1 - minD - 1)):
        for l in range(k + minD + 1, len(RNA) - minI  - 1):
            for j in range(l + minI + 1, min(l + maxI + 1, len(RNA) + 1)):
                if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                    listI += f'I({i},{k},{l},{j})\n'

# create B variabels
for i in range(1, len(RNA) - 1):
    for k in range(1, len(RNA) - 1):
        if k == i+1:
            for l in range(k + minD+1, len(RNA) - 2):
                for j in range(l + 2, min(l + maxB, len(RNA) - 1)):
                    if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                        listB += f'B({i},{k},{l},{j})\n'

        elif k == i-1:
            for l in range(k-minD-1,4,-1):
                for j in range(l - 2, max(l - maxB, 1), -1):
                    if RNA[j-1] + RNA[i-1] in cbp_list and RNA[l-1] + RNA[k-1] in cbp_list:
                        listB += f'B({j},{l},{k},{i})\n'


# create single nucleotide B variabels
# for i in range(1, len(RNA) - 1):
#     for k in range(1, len(RNA) - 1):
#         if k == i+1:
#             for l in range(k + minD, len(RNA) - 2):
#                 for j in range(l + 2, min(l + 2 + B, len(RNA) - 1)):
#                     if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
#                         listB += f'B({i},{k},{l},{j})\n'

#         elif k == i-1:
#             for l in range(k-minD,4,-1):
#                 for j in range(l-2, max(l - 2 - B, 1), -1):
#                     if RNA[j-1] + RNA[i-1] in cbp_list and RNA[l-1] + RNA[k-1] in cbp_list:
#                         listB+= f'B({j},{l},{k},{i})\n'


# print(listB)
# print(listP)
# print(listQ)
# print(listH)
# print(listI)

# print(RNA)

# print(len(listQ.split(', ')))
# print(len(listH.split(', ')))
# print(len(listI.split(', ')))
















# OUT.write("\n")
# OUT.write ("general\n")
# OUT.write("E\n")
# OUT.write  ("end \n")
        
# OUT.close()
