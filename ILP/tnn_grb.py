from binary_variables_grb import *
from stem_parameters import * 
from hairpin_parameters import *
from internal_parameters import *
from initiation_parameters import *
from terminal_mismatch import *
import ilp_parameters_grb

#TODO: stem, hairpin, inernal, bulge classes
#TODO: add multibranch loops, dangling ends and other structural motifs

def num_to_range(num, inMin, inMax, outMin, outMax):
  return outMin + (float(num - inMin) / float(inMax - inMin) * (outMax - outMin))

sc = 1
hc = 1
ic = hc
init = 1
M = 10000
cc = 100

def G_stem(RNA,i,j):
    G = wcf_df.loc[RNA[i-1] + RNA[j-1], RNA[i] + RNA[j-2]]
    return round(G)

def G_F(RNA,i,j):
    if RNA[i-1] + RNA[j-1] == 'AU' or RNA[i-1] + RNA[j-1] == 'UA' or RNA[i-1] + RNA[j-1] == 'GU' or RNA[i-1] + RNA[j-1] == 'UG':        
        G = wcf_initiation*init + wcf_AU_end_penalty
    else: 
        G = wcf_initiation*init
    return round(G)

def G_L(RNA,i,j):
    if RNA[i] + RNA[j-2] == 'AU' or RNA[i] + RNA[j-2] == 'UA' or RNA[i] + RNA[j-2] == 'GU' or RNA[i] + RNA[j-2] == 'UG':
        G = wcf_AU_end_penalty
    else:
        G = 0.0
        # G = G + mismatch_df.loc[RNA[i] + RNA[j-2],RNA[i+1]][RNA[j-3]]
    return round(G)

def G_hairpin(RNA,i,j):

    if j-i-1 == 3:
        G = initiation_df.loc[j-i-1,"hairpin"]
    else:
   
        if j-i-1 <= noH:
            G1 = initiation_df.loc[j-i-1,"hairpin"] + mismatch_df.loc[RNA[i-1] + RNA[j-1],RNA[i]][RNA[j-2]] + spec_GU_clos
            G2 = initiation_df.loc[j-i-1,"hairpin"] + mismatch_df.loc[RNA[i-1] + RNA[j-1],RNA[i]][RNA[j-2]]
        else:
            G1 = M
            G2 = M

        G = G1 if (RNA[i-1] + RNA[j-1] == 'GU' and RNA[i-3] + RNA[i-2] == 'GG') else G2  

        if RNA[i] + RNA[j-2] == 'GA' or RNA[i] + RNA[j-2] == 'UU':
            G = G + hp_mismatch['GA']

        if RNA[i] + RNA[j-2] == 'GG':
            G = G + hp_mismatch['GG']

    return round(G)

def int11(i,k,l,j):
    if k-i-1 == 1 and j-l-1 == 1:
        return True

def G_internal_11(RNA,i,k,l,j):
    return int11_df.loc[RNA[i-1] + RNA[j-1], RNA[i]][RNA[k-1] + RNA[l-1],RNA[l]]

def int12(i,k,l,j):
    if k-i==2 and j-l==3:
        return True
    
def int21(i,k,l,j):
    if k-i==3 and j-l==2:
        return True

def G_internal_12(RNA,i,k,l,j):
    return int12_df.loc[RNA[i-1] + RNA[j-1], RNA[i]][RNA[l],RNA[k-1] + RNA[l-1],RNA[l+1]]

def G_internal_21(RNA,i,k,l,j):
    return int12_df.loc[RNA[l-1] + RNA[k-1], RNA[l]][RNA[i],RNA[j-1] + RNA[i-1],RNA[i+1]]

def int22(i,k,l,j):
    if k-i==3 and j-l==3:
        return True
    
def G_internal_22(RNA,i,k,l,j):
    return int22_df.loc[RNA[i-1] + RNA[j-1],RNA[i]+RNA[l+1]][RNA[k-1] + RNA[l-1],RNA[i+1]+RNA[l]]

def int1n(i,k,l,j):
    if k-i==2 or j-l==2:
        return True
    
def int23(i,k,l,j):
    if (k-i==3 and j-l==4) or (k-i==4 and j-l==3):
        return True
    
def penalty2(RNA,i,k,l,j):
    AU_closure_1 = RNA[i-1] + RNA[j-1] == 'AU' or RNA[i-1] + RNA[j-1] == 'UA'
    GU_closure_1 = RNA[i-1] + RNA[j-1] == 'GU' or RNA[i-1] + RNA[j-1] == 'UG'
    AU_closure_2 = RNA[k-1] + RNA[l-1] == 'AU' or RNA[k-1] + RNA[l-1] == 'UA'  
    GU_closure_2 = RNA[k-1] + RNA[l-1] == 'GU' or RNA[k-1] + RNA[l-1] == 'UG'
    return (AU_closure_1 or GU_closure_1) and (AU_closure_2 or GU_closure_2)

def penalty1(RNA,i,k,l,j):
    AU_closure_1 = RNA[i-1] + RNA[j-1] == 'AU' or RNA[i-1] + RNA[j-1] == 'UA'
    GU_closure_1 = RNA[i-1] + RNA[j-1] == 'GU' or RNA[i-1] + RNA[j-1] == 'UG'
    AU_closure_2 = RNA[k-1] + RNA[l-1] == 'AU' or RNA[k-1] + RNA[l-1] == 'UA'  
    GU_closure_2 = RNA[k-1] + RNA[l-1] == 'GU' or RNA[k-1] + RNA[l-1] == 'UG'
    return (AU_closure_1 or GU_closure_1) or (AU_closure_2 or GU_closure_2)


def G_internal(RNA,i,k,l,j):

    if k-i-1+j-l-1 <= noI:
        common_term = initiation_df.loc[k-i-1+j-l-1,"internal"] + asymmetry * abs(k-i-1-(j-l-1))
    else:
        common_term = M

    if int11(i,k,l,j):
        G = G_internal_11(RNA,i,k,l,j)
        return round(G)
    elif int12(i,k,l,j):
        G = G_internal_12(RNA,i,k,l,j)        
        return round(G)
    elif int21(i,k,l,j):
        G = G_internal_21(RNA,i,k,l,j)        
        return round(G)
    elif int22(i,k,l,j):
        G = G_internal_22(RNA,i,k,l,j)
        return round(G)  
    elif int23(i,k,l,j):        
        if penalty2(RNA,i,k,l,j):
            G = common_term + int23_df.loc[RNA[i-1] + RNA[j-1]][RNA[i] + RNA[j-2]] + int23_df.loc[RNA[l-1] + RNA[k-1]][RNA[l] + RNA[k-2]] + 2*AU_end_penalty
            return round(G) 
        elif penalty1(RNA,i,k,l,j):
            G = common_term + int23_df.loc[RNA[i-1] + RNA[j-1]][RNA[i] + RNA[j-2]] + int23_df.loc[RNA[l-1] + RNA[k-1]][RNA[l] + RNA[k-2]] + AU_end_penalty
            return round(G) 
        else:
            G = common_term + int23_df.loc[RNA[i-1] + RNA[j-1]][RNA[i] + RNA[j-2]] + int23_df.loc[RNA[l-1] + RNA[k-1]][RNA[l] + RNA[k-2]]
            return round(G)   
    elif int1n(i,k,l,j):
        if penalty2(RNA,i,k,l,j):
            G = common_term + 2*AU_end_penalty
            return round(G) 
        elif penalty1(RNA,i,k,l,j):
            G = common_term + AU_end_penalty         
            return round(G) 
        else:
            G = common_term
            return round(G)
    else:
        if penalty2(RNA,i,k,l,j):
            G = common_term + intnn_df.loc[RNA[i-1] + RNA[j-1]][RNA[i] + RNA[j-2]] + intnn_df.loc[RNA[l-1] + RNA[k-1]][RNA[l] + RNA[k-2]] + 2*AU_end_penalty
            return round(G) 
        elif penalty1(RNA,i,k,l,j):
            G = common_term + intnn_df.loc[RNA[i-1] + RNA[j-1]][RNA[i] + RNA[j-2]] + intnn_df.loc[RNA[l-1] + RNA[k-1]][RNA[l] + RNA[k-2]] + AU_end_penalty
            return round(G) 
        else:
            G = common_term + intnn_df.loc[RNA[i-1] + RNA[j-1]][RNA[i] + RNA[j-2]] + intnn_df.loc[RNA[l-1] + RNA[k-1]][RNA[l] + RNA[k-2]]
            return round(G)

def G_bulge(RNA,i,k,l,j):
    if k-i-1+j-l-1 <= noB:
        if k==i+1:
            if j-l-1 == 1:
                G = initiation_df.loc[j-l-1,"bulge"] + wcf_df.loc[RNA[i-1] + RNA[j-1], RNA[k-1] + RNA[l-1]] + Cbulge*(RNA[l] == "C" and (RNA[l-1] == "C" or RNA[j-1] == "C")) - RT*np.log(3)
            if j-l-1 > 1 and j-l-1 <= 6:
                G = initiation_df.loc[j-l-1,"bulge"]
            if j-l-1 > 6:
                G = initiation_df.loc[j-l-1,"bulge"] + 1.75*RT*np.log((j-l-1)/6)
        elif j==l+1:
            if k-i-1 == 1: 
                G = initiation_df.loc[k-i-1,"bulge"] + wcf_df.loc[RNA[i-1] + RNA[j-1], RNA[k-1] + RNA[l-1]] + Cbulge*(RNA[k] == "C" and (RNA[k-1] == "C" or RNA[i-1] == "C")) - RT*np.log(3)
            if k-i-1 > 1 and k-i-1 <= 6:
                G = initiation_df.loc[k-i-1,"bulge"]
            if k-i-1 > 6:
                G = initiation_df.loc[k-i-1,"bulge"] + 1.75*RT*np.log((k-i-1)/6)
    else:
        G = M
    return round(G)

def G_multi(i,i1,j1,i2,j2,j):
    if (i1-i-1 > maxM) or (i2-j1-1 > maxM) or (j2-j-1 > maxM):
        G = M
    else:
        G = cc*(ilp_parameters_grb.c*(i1-i-1+i2-j1-1+j-j2-1) + b*2)
    return round(G)

#RNA = 'AACCAUGUCAGGUCCGGAAGGAAGCAGCAU'
# RNA = 'CAGACGCGGAGUG'
# i=2
# k=5
# l=8
# j=12
# print(G_internal(RNA,i,k,l,j))


# common_term = initiation_df.loc[k-i-1+j-l-1,"internal"] + asymmetry * abs(k-i-1-(j-l-1))

# print(initiation_df.loc[k-i-1+j-l-1,"internal"])
# G = common_term + int23_df.loc[RNA[i-1] + RNA[j-1]][RNA[i] + RNA[j-2]] + int23_df.loc[RNA[l-1] + RNA[k-1]][RNA[l] + RNA[k-2]] + AU_end_penalty

# print(initiation_df.loc[k-i-1+j-l-1,"internal"])
# print(asymmetry * abs(k-i-1-(j-l-1)))
# print(int23_df.loc[RNA[i-1] + RNA[j-1]][RNA[i] + RNA[j-2]])
# print(int23_df.loc[RNA[l-1] + RNA[k-1]][RNA[l] + RNA[k-2]])
# print(AU_end_penalty)







