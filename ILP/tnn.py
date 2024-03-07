from binary_variables import *
from stem_parameters import * 
from hairpin_parameters import *
from internal_parameters import *
from initiation_parameters import *
from terminal_mismatch import *

#TODO: stem, hairpin, inernal, bulge classes
#TODO: add multibranch loops, dangling ends and other structural motifs

def num_to_range(num, inMin, inMax, outMin, outMax):
  return outMin + (float(num - inMin) / float(inMax - inMin) * (outMax - outMin))

sc = 1
hc = 1
ic = hc
init = 1

def G_stem(i,j):
    G = f'{wcf_df.loc[RNA[i-1] + RNA[j-1], RNA[i] + RNA[j-2]]}'
    return G

def G_F(i,j):
    if RNA[i-1] + RNA[j-1] == 'AU' or RNA[i-1] + RNA[j-1] == 'UA' or RNA[i-1] + RNA[j-1] == 'GU' or RNA[i-1] + RNA[j-1] == 'UG':        
        G = f'{(wcf_initiation*init + wcf_AU_end_penalty)}' 
    else: 
        G = f'{wcf_initiation*init}'
    return G

def G_L(i,j):
    if RNA[i] + RNA[j-2] == 'AU' or RNA[i] + RNA[j-2] == 'UA' or RNA[i-1] + RNA[j-1] == 'GU' or RNA[i-1] + RNA[j-1] == 'UG':
        G = f'{wcf_AU_end_penalty}' 
    else:
        G = f'{0.0}'

    return G

def G_hairpin(i,j):

    G1 = initiation_df.loc[j-i-1,"hairpin"] + mismatch_df.loc[RNA[i-1] + RNA[j-1],RNA[i]][RNA[j-2]] + spec_GU_clos
    G2 = initiation_df.loc[j-i-1,"hairpin"] + mismatch_df.loc[RNA[i-1] + RNA[j-1],RNA[i]][RNA[j-2]]

    G = f'{G1}' if RNA[i-1] + RNA[j-1] == 'GU' else f'{G2}'     
    return G

def int11(i,k,l,j):
    if k-i-1 == 1 and j-l-1 == 1:
        return True

def G_internal_11(i,k,l,j):
    # print(f'{RNA[i-1] + RNA[j-1]},{RNA[i]}::{RNA[k-1] + RNA[l-1]},{RNA[l]}')
    return int11_df.loc[RNA[i-1] + RNA[j-1], RNA[i]][RNA[k-1] + RNA[l-1],RNA[l]]

def int12(i,k,l,j):
    if k-i-1==1 and j-l==3:
        return True

def G_internal_12(i,k,l,j):
    return int12_df.loc[RNA[i-1] + RNA[j-1], RNA[i]][RNA[l],RNA[k-1] + RNA[l-1],RNA[l+1]]

def int22(i,k,l,j):
    if k-i==3 and j-l==3:
        return True
    
def G_internal_22(i,k,l,j):
    return int22_df.loc[RNA[i-1] + RNA[j-1],RNA[i]+RNA[l+1]][RNA[k-1] + RNA[l-1],RNA[i+1]+RNA[l]]

def G_internal(i,k,l,j):

    if int11(i,k,l,j):
        # print("1x1")
        return G_internal_11(i,k,l,j)
    elif int12(i,k,l,j):
        # print("1x2")
        return G_internal_12(i,k,l,j)
    elif int22(i,k,l,j):
        # print("2x2")
        return G_internal_22(i,k,l,j)
    else:
        # print("n1xn2")
        if (RNA[i-1] + RNA[j-1] == 'GU' or RNA[i-1] + RNA[j-1] == 'AU' and RNA[k-1] + RNA[l-1] == 'GU' or RNA[k-1] + RNA[l-1] == 'AU'):
            return f'{initiation_df.loc[k-i-1+j-l-1,"internal"] + asymmetry * abs(k-i-1-(j-l-1)) + mismatch_df.loc[RNA[i-1] + RNA[j-1],RNA[i]][RNA[j-2]] + mismatch_df.loc[RNA[k-1] + RNA[l-1],RNA[k-2]][RNA[l]] + 2*AU_end_penalty}' 
        elif (RNA[i-1] + RNA[j-1] == 'GU' or RNA[i-1] + RNA[j-1] == 'AU' or RNA[k-1] + RNA[l-1] == 'GU' or RNA[k-1] + RNA[l-1] == 'AU'):
            return f'{initiation_df.loc[k-i-1+j-l-1,"internal"] + asymmetry * abs(k-i-1-(j-l-1)) + mismatch_df.loc[RNA[i-1] + RNA[j-1],RNA[i]][RNA[j-2]] + mismatch_df.loc[RNA[k-1] + RNA[l-1],RNA[k-2]][RNA[l]] + AU_end_penalty}' 
        else:
            return f'{initiation_df.loc[k-i-1+j-l-1,"internal"] + asymmetry * abs(k-i-1-(j-l-1)) + mismatch_df.loc[RNA[i-1] + RNA[j-1],RNA[i]][RNA[j-2]] + mismatch_df.loc[RNA[k-1] + RNA[l-1],RNA[k-2]][RNA[l]]}' 

def G_bulge(i,k,l,j):

    if k==i+1:
        if j-l-1 == 1:
            G = f'{initiation_df.loc[j-l-1,"bulge"] + wcf_df.loc[RNA[i-1] + RNA[j-1], RNA[k-1] + RNA[l-1]] + Cbulge*(RNA[l] == "C") - RT*np.log(3)}'
        if j-l-1 > 1 and j-l-1 <= 6:
            G = f'{initiation_df.loc[j-l-1,"bulge"]}'
        if j-l-1 > 6:
            G = f'{initiation_df.loc[j-l-1,"bulge"] + 1.75*RT*np.log((j-l-1)/6)}'
    elif j==l+1:
        if k-i-1 == 1: 
            G = f'{initiation_df.loc[k-i-1,"bulge"] + wcf_df.loc[RNA[i-1] + RNA[j-1], RNA[k-1] + RNA[l-1]] + Cbulge*(RNA[k] == "C") - RT*np.log(3)}'

        if k-i-1 > 1 and k-i-1 <= 6:
            G = f'{initiation_df.loc[k-i-1,"bulge"]}'
        if k-i-1 > 6:
            G = f'{initiation_df.loc[k-i-1,"bulge"] + 1.75*RT*np.log((k-i-1)/6)}'

    return G
    

# print(G_stem(1,9))
# print(G_hairpin(1,5))
# print(G_internal(1,3,17,20))
# print(G_bulge(3,4,8,10))
# print(G_bulge(11,12,19,21))
# print(G_bulge(2,3,18,22))

# print(G_hairpin(13,17))
# print(G_internal(5,8,23,25))
# print(G_internal(10,12,18,21))

print(G_hairpin(11,16))
print(G_internal(9, 11, 16, 18))
print(G_stem(8,19))

print(G_internal(5,8,19,21))

print(G_stem(4,22))
print(G_stem(3,23))
print(G_stem(2,24))
print(G_stem(1,25))


