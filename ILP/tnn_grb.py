from binary_variables_grb import *
from stem_parameters import * 
from hairpin_parameters import *
from internal_parameters import *
from initiation_parameters import *
from terminal_mismatch import *

from icecream import ic

#TODO: stem, hairpin, inernal, bulge classes
#TODO: add multibranch loops, dangling ends and other structural motifs

def num_to_range(num, inMin, inMax, outMin, outMax):
  return outMin + (float(num - inMin) / float(inMax - inMin) * (outMax - outMin))

sc = 1
hc = 1
ic = hc
init = 1
M = 10000

def G_stem(i,j):
    G = wcf_df.loc[RNA[i-1] + RNA[j-1], RNA[i] + RNA[j-2]]
    return round(G)

def G_F(i,j):
    if RNA[i-1] + RNA[j-1] == 'AU' or RNA[i-1] + RNA[j-1] == 'UA' or RNA[i-1] + RNA[j-1] == 'GU' or RNA[i-1] + RNA[j-1] == 'UG':        
        G = wcf_initiation*init + wcf_AU_end_penalty
    else: 
        G = wcf_initiation*init
    return round(G)

def G_L(i,j):
    if RNA[i] + RNA[j-2] == 'AU' or RNA[i] + RNA[j-2] == 'UA' or RNA[i] + RNA[j-2] == 'GU' or RNA[i] + RNA[j-2] == 'UG':
        G = wcf_AU_end_penalty
    else:
        G = 0.0

    # G = G + mismatch_df.loc[RNA[i] + RNA[j-2],RNA[i+1]][RNA[j-3]]

    return round(G)

def G_hairpin(i,j):

    # G1 = initiation_df.loc[j-i-1,"hairpin"] + mismatch_df.loc[RNA[i-1] + RNA[j-1],RNA[i]][RNA[j-2]] + spec_GU_clos
    # G2 = initiation_df.loc[j-i-1,"hairpin"] + mismatch_df.loc[RNA[i-1] + RNA[j-1],RNA[i]][RNA[j-2]]
   
    if j-i-1 <= noH:
        G1 = initiation_df.loc[j-i-1,"hairpin"] + mismatch_df.loc[RNA[i-1] + RNA[j-1],RNA[i]][RNA[j-2]] + spec_GU_clos
        G2 = initiation_df.loc[j-i-1,"hairpin"] + mismatch_df.loc[RNA[i-1] + RNA[j-1],RNA[i]][RNA[j-2]]
    else:
        G1 = M + spec_GU_clos
        G2 = M

    G = G1 if RNA[i-1] + RNA[j-1] == 'GU' else G2  

    if RNA[i] + RNA[j-2] == 'GA' or RNA[i] + RNA[j-2] == 'UU':
        G = G + hp_mismatch['GA']

    if RNA[i] + RNA[j-2] == 'GG':
        G = G + hp_mismatch['GG']

    return round(G)

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

def int1n(i,k,l,j):
    if k-i==2 or j-l==2:
        return True
    
def int23(i,k,l,j):
    if k-i==3 and j-l==4:
        return True
    
def penalty2(i,k,l,j):
    return (RNA[i-1] + RNA[j-1] == 'GU' or RNA[i-1] + RNA[j-1] == 'AU' and RNA[k-1] + RNA[l-1] == 'GU' or RNA[k-1] + RNA[l-1] == 'AU')

def penalty1(i,k,l,j):
    return (RNA[i-1] + RNA[j-1] == 'GU' or RNA[i-1] + RNA[j-1] == 'AU' or RNA[k-1] + RNA[l-1] == 'GU' or RNA[k-1] + RNA[l-1] == 'AU')


def G_internal(i,k,l,j):

    if k-i-1+j-l-1 <= noI:
        common_term = initiation_df.loc[k-i-1+j-l-1,"internal"] + asymmetry * abs(k-i-1-(j-l-1))
    else:
        common_term = M + asymmetry * abs(k-i-1-(j-l-1))

    if int11(i,k,l,j):
        # print("1x1")
        G = G_internal_11(i,k,l,j)
        return round(G)
    elif int12(i,k,l,j):
        # print("1x2")
        G = G_internal_12(i,k,l,j)
        return round(G)
    elif int22(i,k,l,j):
        # print("2x2")
        G = G_internal_22(i,k,l,j)
        return round(G)  
    elif int23(i,k,l,j):        
        if penalty2(i,k,l,j):
            G = common_term + int23_df.loc[RNA[i-1] + RNA[j-1]][RNA[i] + RNA[j-2]] + int23_df.loc[RNA[l-1] + RNA[k-1]][RNA[l] + RNA[k-2]] + 2*AU_end_penalty
            return round(G) 
        elif penalty1(i,k,l,j):
            G = common_term + int23_df.loc[RNA[i-1] + RNA[j-1]][RNA[i] + RNA[j-2]] + int23_df.loc[RNA[l-1] + RNA[k-1]][RNA[l] + RNA[k-2]] + AU_end_penalty
            return round(G) 
        else:
            G = common_term + int23_df.loc[RNA[i-1] + RNA[j-1]][RNA[i] + RNA[j-2]] + int23_df.loc[RNA[l-1] + RNA[k-1]][RNA[l] + RNA[k-2]]
            return round(G)   
    elif int1n(i,k,l,j):
        # print("1xn")
        if penalty2(i,k,l,j):
            G = common_term + 2*AU_end_penalty
            return round(G) 
        elif penalty1(i,k,l,j):
            G = common_term + AU_end_penalty         
            return round(G) 
        else:
            G = common_term
            return round(G)
    else:
        # print("n1xn2")
        if penalty2(i,k,l,j):
            G = common_term + intnn_df.loc[RNA[i-1] + RNA[j-1]][RNA[i] + RNA[j-2]] + intnn_df.loc[RNA[l-1] + RNA[k-1]][RNA[l] + RNA[k-2]] + 2*AU_end_penalty
            return round(G) 
        elif penalty1(i,k,l,j):
            G = common_term + intnn_df.loc[RNA[i-1] + RNA[j-1]][RNA[i] + RNA[j-2]] + intnn_df.loc[RNA[l-1] + RNA[k-1]][RNA[l] + RNA[k-2]] + AU_end_penalty
            return round(G) 
        else:
            G = common_term + intnn_df.loc[RNA[i-1] + RNA[j-1]][RNA[i] + RNA[j-2]] + intnn_df.loc[RNA[l-1] + RNA[k-1]][RNA[l] + RNA[k-2]]
            return round(G)

def G_bulge(i,k,l,j):

    if k-i-1+j-l-1 <= noB:

        if k==i+1:
            if j-l-1 == 1:
                G = initiation_df.loc[j-l-1,"bulge"] + wcf_df.loc[RNA[i-1] + RNA[j-1], RNA[k-1] + RNA[l-1]] + Cbulge*(RNA[l] == "C") - RT*np.log(3)
            if j-l-1 > 1 and j-l-1 <= 6:
                G = initiation_df.loc[j-l-1,"bulge"]
            if j-l-1 > 6:
                G = initiation_df.loc[j-l-1,"bulge"] + 1.75*RT*np.log((j-l-1)/6)
        elif j==l+1:
            if k-i-1 == 1: 
                G = initiation_df.loc[k-i-1,"bulge"] + wcf_df.loc[RNA[i-1] + RNA[j-1], RNA[k-1] + RNA[l-1]] + Cbulge*(RNA[k] == "C") - RT*np.log(3)
            if k-i-1 > 1 and k-i-1 <= 6:
                G = initiation_df.loc[k-i-1,"bulge"]
            if k-i-1 > 6:
                G = initiation_df.loc[k-i-1,"bulge"] + 1.75*RT*np.log((k-i-1)/6)
    else:
        G = M

    return round(G)

# 1Q9A
# generated -634
# (((((..((........)).)))))
# print(G_stem(1, 25))
# print(G_stem(2, 24))
# print(G_stem(3, 23))
# print(G_stem(4, 22))
# print(G_stem(8, 19))
# print(G_hairpin(9, 18))
# print(G_internal(5,8,19,21))
# print('\\\TOTAL///')
# print(G_stem(1, 25)+G_stem(2, 24)+G_stem(3, 23)+G_stem(4, 22)+G_stem(8, 19)+G_hairpin(9, 18)+G_internal(5,8,19,21))


# alternative -594 :: set H(11,16) == 1
# (((((..((.(....).)).)))))
# print(G_stem(1, 25))
# print(G_stem(2, 24))
# print(G_stem(3, 23))
# print(G_stem(4, 22))
# print(G_stem(8, 19))
# print(G_hairpin(11, 16))
# print(G_internal(5,8,19,21))
# print(G_internal(9,11,16,18))
# print('\\\TOTAL///')
# print(G_stem(1, 25)+G_stem(2, 24)+G_stem(3, 23)+G_stem(4, 22)+G_stem(8, 19)+G_hairpin(11, 16)+G_internal(5,8,19,21)+G_internal(9,11,16,18))


# # 3sn2 test

# # reference MFE
# print(G_stem(1,29))
# print(G_stem(2,28))
# print(G_stem(3,27))
# print(G_stem(4,26))
# print(G_stem(5,25))
# print(G_bulge(6,8,23,24))
# print(G_stem(8,23))
# print(G_stem(9,22))
# print(G_stem(10,21))
# print(G_stem(11,20))
# print(G_bulge(12,13,17,19))
# print(G_hairpin(13,17))

# print('TOTAL REFERENCE')
# print(f'{float(G_stem(1,29)) + float(G_stem(2,28)) + float(G_stem(3,27)) + float(G_stem(4,26)) + float(G_stem(5,25)) + float(G_bulge(6,8,23,24)) + float(G_stem(8,23)) + float(G_stem(9,22)) + float(G_stem(10,21))+float(G_stem(11,20)) + float(G_bulge(12,13,17,19)) + float(G_hairpin(13,17))}')

# # generated MFE
# print(G_stem(1,29))
# print(G_stem(2,28))
# print(G_stem(3,27))
# print(G_stem(4,26))
# print(G_internal(5,8,23,25))
# print(G_stem(8,23))
# print(G_stem(9,22))
# print(G_stem(10,21))
# print(G_stem(11,20))
# print(G_hairpin(12,19))

# print('TOTAL GENERATED')
# print(f'{float(G_stem(1,29)) + float(G_stem(2,28)) + float(G_stem(3,27)) + float(G_stem(4,26)) + float(G_stem(8,23)) + float(G_stem(9,22)) + float(G_stem(10,21))+float(G_stem(11,20)) + float(G_internal(5,8,23,25)) + float(G_hairpin(12,19))}')




