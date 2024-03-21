from binary_variables import *
from tnn import *

def num_to_range(num, inMin, inMax, outMin, outMax):
  return outMin + (float(num - inMin) / float(inMax - inMin) * (outMax
                  - outMin))

def stemParams(RNA):
    g = {}

    for i in range(1, len(RNA)):  
        for j in range(i + minD + 1, len(RNA) + 1): 
            ip1 = i + 1
            jm1 = j - 1
            if legal(i,j):
                if legal(ip1,jm1) and ip1 < jm1:

                    g[f'Q{i}',f'Q{j}'] = G_stem(i,j)
                    
    return g

def flParams(RNA):
    g = {}

    for i in range(1, len(RNA)):  
        for j in range(i + minD + 1, len(RNA) + 1): 
            if legal(i,j):
                if legal(i+1,j-1) and i+1 < j-1:
                    g[f'F{i}',f'F{j}'] = G_F(i,j)
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

def normalize(RNA):
    q_params = {k: float(v) for k, v in stemParams(RNA).items()}
    fl_params = {k: float(v) for k, v in flParams(RNA).items()}
    h_params = {k: float(v) for k, v in hairpinParams(RNA).items()}
    i_params = {k: float(v) for k, v in internalParams(RNA).items()}
    b_params = {k: float(v) for k, v in bulgeParams(RNA).items()}


    params = q_params | fl_params | h_params | i_params | b_params

    # inMax = max(params.values())
    # inMin = min(params.values())
    # outMin = 0.0
    # outMax = 100.0

    # print(inMin)
    # print(inMax)

    # normalized = {k: round(num_to_range(float(v), float(inMin), float(inMax), outMin, outMax)) for k, v in params.items()}

    normalized = {k: round(-float(v)) for k, v in params.items()}

    return normalized

# pretty(normalize(RNA))

def stemTerm(RNA):
    objective = ''
    norm = normalize(RNA)

    for i in range(1, len(RNA)):  
        for j in range(i + minD + 1, len(RNA)+1): 
            ip1 = i + 1
            jm1 = j - 1
            if legal(i,j):
                if legal(ip1,jm1) and ip1 < jm1:                    
                    objective += f" + {norm[f'Q{i}',f'Q{j}']} Q({i},{j})"

    return objective

def fTerm(RNA):
    objective = ''
    norm = normalize(RNA)

    for i in range (1, len(RNA)): 
        for j in range(i + minD + 1,  len(RNA)+1): 
            if legal(i,j) and legal(i+1,j-1):  
                objective += f" + {norm[f'F{i}',f'F{j}']} F({i},{j})" 

    return objective

def lTerm(RNA):
    objective = ''
    norm = normalize(RNA)

    for i in range (1, len(RNA)): 
        for j in range(i + minD + 1,  len(RNA)+1): 
            if i > 1 and j < len(RNA):
                if legal(i,j) and legal(i+1,j-1):                
                    objective += f" + {norm[f'L{i}',f'L{j}']} L({i},{j})"

    return objective

def hairpinTerm(RNA):
    objective = ''
    norm = normalize(RNA)

    for i in range(1,len(RNA) - minD - 1):
        for j in range(i + minD + 1, min(i + maxH + 1, len(RNA) + 1)):
            if RNA[i-1] + RNA[j-1] in cbp_list:     
                objective += f" + {norm[f'H{i}',f'H{j}']} H({i},{j})"
    return objective

def internalTerm(RNA):
    objective = ''
    norm = normalize(RNA)

    for i in range(1, len(RNA) - minI - 1 - minD - 1 - minI - 1):
        for k in range(i + minI + 1, min(i + maxI + 1, len(RNA) - minI - 1 - minD - 1)):
            for l in range(k + minD + 1, len(RNA) - minI  - 1):
                for j in range(l + minI + 1, min(l + maxI + 1, len(RNA) + 1)):
                    if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:                        
                        objective += f" + {norm[f'I{i}',f'I{k}',f'I{l}',f'I{j}']} I({i},{k},{l},{j})"

    return objective

def bulgeTerm(RNA):
    objective = ''
    norm = normalize(RNA)

    for i in range(1, len(RNA) - 1):
        for k in range(1, len(RNA) - 1):
            if k == i+1:
                for l in range(k + minD + 1, len(RNA) - 2):
                    for j in range(l + 2, min(l + maxB, len(RNA) - 1)):
                        if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
                            objective += f" + {norm[f'B{i}',f'B{k}',f'B{l}',f'B{j}']} B({i},{k},{l},{j})"

            elif k == i-1:
                for l in range(k-minD-1, 4, -1):
                    for j in range(l - 2, max(l - maxB, 1), -1):
                        if RNA[j-1] + RNA[i-1] in cbp_list and RNA[l-1] + RNA[k-1] in cbp_list:
                            objective += f" + {norm[f'B{j}',f'B{l}',f'B{k}',f'B{i}']} B({j},{l},{k},{i})"
    
    return objective

# print(stemTerm(RNA))
# print(fTerm(RNA))
# print(lTerm(RNA))
# print(hairpinTerm(RNA))
# print(internalTerm(RNA))
# print(bulgeTerm(RNA))