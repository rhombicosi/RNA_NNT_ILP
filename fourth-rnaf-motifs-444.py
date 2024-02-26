import sys
import math
            
from score_matrix import bp_df, stack_df 
# bp_df=bp_df.set_index('Desc')
     


with open('6n5n', 'r') as file:
    RNA = file.read().rstrip()    
    print (RNA)

length = len(RNA)
print (length)
    
minD = 3

loop6 = 96997.0
loop6_s = 78324
loop6_pq = 18673
loop6_q = 16547
loop6_p = 2126
    
loop8 = 112139.0
loop8_s = 106865
loop8_pq = 5274.0
loop8_p = 1239
loop8_q = 4035

listP = ""      # listP will sum the edge-to-edge pairs variables
listQ = ""      # listQ will sum the face-to-face stacking variables
listS = ""

OUT = open('6n5n-motifs-loop6.lp', "w")
OUT.write ("Maximize \n")
OUT.write ("E \n")

length_for_range = length + 1
weights = [[0 for x in range(1, length_for_range + 1)] for y in range(1, length_for_range + 1)]

# 6n5n
# loop8
a_start = 7#46
a_end = 12#51

b_start = 22#67
b_end = 25#70

c_start = 106#86
c_end = 109#89

# loop6
# a_start = 46
# a_end = 51

# b_start = 67
# b_end = 70

# c_start = 86
# c_end = 89

A = range(a_start,a_end)
B = range(b_start,b_end)
C = range(c_start,c_end)

RNA_CUT = RNA[a_start-1:a_end] + RNA[b_start-1:b_end] + RNA[c_start-1:c_end]
print(RNA_CUT)
print(len(RNA_CUT))

# RNA = RNA_CUT

# a_start = 1 #46
# a_end = 6 #51


# b_start = 7 #67
# b_end = 10 #70

# c_start = 11 #86
# c_end = 14 #89

# all edge2edge contacts
NP_list = ['AU','UA','CG','GC','GU','UG']
def listp_construct(a_start,a_end,b_start,b_end,listP):
    
    for i in range(a_start+1, a_end):
        for j in range(b_start+1, b_end):
            if (j-i) > minD:
                
                # edge2edge contacts weights
                ww_weight = bp_df.loc['WW_cis', RNA[i-1] + RNA[j-1]]
                hh_weight = bp_df.loc['HH_cis', RNA[i-1] + RNA[j-1]]
                ss_weight = bp_df.loc['SS_cis', RNA[i-1] + RNA[j-1]]
                
                ws_weight = bp_df.loc['WS_cis', RNA[i-1] + RNA[j-1]]
                wh_weight = bp_df.loc['WH_cis', RNA[i-1] + RNA[j-1]]
                hs_weight = bp_df.loc['HS_cis', RNA[i-1] + RNA[j-1]]
                
                if not math.isnan(float(ww_weight)):
                    if RNA[i-1] + RNA[j-1] not in NP_list:                        
                        listP = listP + f' + {ww_weight} WW({i},{j})'
                if not math.isnan(float(hh_weight)):
                    listP = listP + f' + {hh_weight} HH({i},{j})'
                if not math.isnan(float(ss_weight)):
                    listP = listP + f' + {ss_weight} SS({i},{j})'
                    
                if not math.isnan(float(wh_weight)):
                    listP = listP + f' + {wh_weight} WH({i},{j})'
                    listP = listP + f' + {wh_weight} WH({j},{i})'
                if not math.isnan(float(ws_weight)):
                    listP = listP + f' + {ws_weight} WS({i},{j})'
                    listP = listP + f' + {ws_weight} WS({j},{i})'
                if not math.isnan(float(hs_weight)):
                    listP = listP + f' + {hs_weight} HS({i},{j})'
                    listP = listP + f' + {hs_weight} HS({j},{i})'                            
                        
    return listP
                            
listP = listp_construct(a_start,a_end,b_start,b_end,listP)
listP = listp_construct(b_start,b_end,c_start,c_end,listP)
listP = listp_construct(c_start,c_end,a_start,a_end,listP)

# print(listP)

# face2face contacts btwn 2 different segments
def listq_construct(a_start,a_end,b_start,b_end,listQ):
    
    for i in range(a_start+1, a_end):
        for j in range(b_start+1, b_end):
            if (j-i) > minD:
                                        
                # edge-to-edge contacts weights
                fk_weight = stack_df.loc['>>', RNA[i-1] + RNA[j-1]]
                kf_weight = stack_df.loc['<<', RNA[i-1] + RNA[j-1]]
                ff_weight = stack_df.loc['><', RNA[i-1] + RNA[j-1]]
                kk_weight = stack_df.loc['<>', RNA[i-1] + RNA[j-1]]
                
                
                if not math.isnan(float(fk_weight)):
                    listQ = listQ + f' + {fk_weight} FK({i},{j})'
                if not math.isnan(float(kf_weight)):
                    listQ = listQ + f' + {kf_weight} KF({i},{j})'
                if not math.isnan(float(ff_weight)):
                    listQ = listQ + f' + {ff_weight} FF({i},{j})'
                if not math.isnan(float(kk_weight)):
                    listQ = listQ + f' + {kk_weight} KK({i},{j})'
                        
    return listQ          
                    
listQ = listq_construct(a_start,a_end,b_start,b_end,listQ)
listQ = listq_construct(b_start,b_end,c_start,c_end,listQ)
listQ = listq_construct(c_start,c_end,a_start,a_end,listQ)

# print(listQ)

# face2face contacts allowed within the same segment
def lists_construct(a_start,a_end, listS):
    
    for i in range(a_start, a_end):
        fk_weight = stack_df.loc['>>', RNA[i] + RNA[i+1]]
        
        if not math.isnan(float(fk_weight)):
            listS = listS + f' + {fk_weight} FK({i},{i+1})'
                        
    return listS

listS = lists_construct(a_start,a_end, listS)
listS = lists_construct(b_start,b_end, listS)
listS = lists_construct(c_start,c_end, listS)

OUT.write ("\nsuch that \n")
OUT.write (listP + listQ + listS + " - E = 0 \n\n")

# The following lines generate the inequalities to ensure
# that each edge is paired to at most one other edge

# W edge inequalities
def w_inequalities(a_start,a_end,b_start,b_end,c_start,c_end):
    for i in range(a_start + 1, a_end):
        inequality = ""    
        for j in range(b_start + 1, b_end):
            if RNA[i-1] + RNA[j-1] not in NP_list:
                inequality = inequality + f' + WW({i},{j})'
            inequality = inequality + f' + WH({i},{j})'
            inequality = inequality + f' + WS({i},{j})'
            
            # inequality = inequality + f' + WH({j},{i})'
            # inequality = inequality + f' + WS({j},{i})'
            
        for j in range(c_start + 1, c_end):
            if RNA[j-1] + RNA[i-1] not in NP_list:
                inequality = inequality + f' + WW({j},{i})'
            inequality = inequality + f' + WH({i},{j})'
            inequality = inequality + f' + WS({i},{j})'
            
            # inequality = inequality + f' + WH({j},{i})'
            # inequality = inequality + f' + WS({j},{i})'
        
        inequality = inequality + ' <= 1'
        OUT.write (inequality)           		
        OUT.write ("\n\n")

w_inequalities(a_start,a_end,b_start,b_end,c_start,c_end)
w_inequalities(b_start,b_end,c_start,c_end,a_start,a_end)
w_inequalities(c_start,c_end,a_start,a_end,b_start,b_end)
        
# H edge inequalities
def  h_inequalities(a_start,a_end,b_start,b_end,c_start,c_end):
    for i in range(a_start + 1, a_end):
        inequality = ""    
        for j in range(b_start + 1, b_end):        
      
            inequality = inequality + f' + HH({i},{j})'
            inequality = inequality + f' + WH({j},{i})'
            inequality = inequality + f' + HS({i},{j})'
        
        for j in range(c_start + 1, c_end):        
      
            inequality = inequality + f' + HH({j},{i})'
            inequality = inequality + f' + WH({j},{i})'
            inequality = inequality + f' + HS({i},{j})'
                        
        inequality = inequality + ' <= 1'
        OUT.write (inequality)           		
        OUT.write ("\n\n")
            
h_inequalities(a_start,a_end,b_start,b_end,c_start,c_end)
h_inequalities(b_start,b_end,c_start,c_end,a_start,a_end)
h_inequalities(c_start,c_end,a_start,a_end,b_start,b_end)

# S edge inequalities
def s_inequalities(a_start,a_end,b_start,b_end,c_start,c_end):
    for i in range(a_start + 1, a_end):
        inequality = ""    
        for j in range(b_start + 1, b_end):
            
            inequality = inequality + f' + SS({i},{j})'        
            inequality = inequality + f' + WS({j},{i})'
            inequality = inequality + f' + HS({j},{i})'
            
            # inequality = inequality + f' + WS({j},{i})'
            # inequality = inequality + f' + HS({j},{i})'
            
        for j in range(c_start + 1,c_end):
            
            inequality = inequality + f' + SS({j},{i})'        
            inequality = inequality + f' + WS({j},{i})'
            inequality = inequality + f' + HS({j},{i})'
            
            # inequality = inequality + f' + WS({j},{i})'
            # inequality = inequality + f' + HS({j},{i})'
                        
        inequality = inequality + ' <= 1'
        OUT.write (inequality)           		
        OUT.write ("\n\n")
            
        # print(f'inequality #{i}: {inequality}')
            
s_inequalities(a_start,a_end,b_start,b_end,c_start,c_end)
s_inequalities(b_start,b_end,c_start,c_end,a_start,a_end)
s_inequalities(c_start,c_end,a_start,a_end,b_start,b_end)


# stackings inequalities for 3' (F) face btwn two different segments
# stackings inequalities to ensure one stacking connection btwn two nucleotides
def n_inequalities(a_start,a_end,b_start,b_end,c_start,c_end):
    for i in range(a_start + 1, a_end):  
            inequality = ""        
            for j in range(b_start + 1, b_end):          
                inequality = inequality + f' + FF({i},{j})'
                inequality = inequality + f' + FK({i},{j})'  
                inequality = inequality + f' + KF({i},{j})'
                inequality = inequality + f' + KK({i},{j})'
                    
            for j in range(c_start + 1, c_end):          
                inequality = inequality + f' + FF({j},{i})'
                inequality = inequality + f' + KF({j},{i})'  
                inequality = inequality + f' + FK({j},{i})'
                inequality = inequality + f' + KK({j},{i})'
                
            inequality = inequality + ' <= 1'
            OUT.write (inequality)
            OUT.write ("\n\n")

n_inequalities(a_start,a_end,b_start,b_end,c_start,c_end)
n_inequalities(b_start,b_end,c_start,c_end,a_start,a_end)
n_inequalities(c_start,c_end,a_start,a_end,b_start,b_end)

# stackings inequalities for max one stacking for f(3') face 
def f_inequalities(a_start,a_end,b_start,b_end,c_start,c_end):
    for i in range(a_start + 1, a_end):  
        inequality = ""        
        inequality = inequality + f' + FK({i},{i+1})'
        
        for j in range(b_start + 1, b_end):  
            inequality = inequality + f' + FK({i},{j})'
            inequality = inequality + f' + FF({i},{j})'
            
        for j in range(c_start + 1, c_end):  
            inequality = inequality + f' + KF({j},{i})'
            inequality = inequality + f' + FF({j},{i})'
        
        inequality = inequality + ' <= 1\n'
        OUT.write (inequality)
        
f_inequalities(a_start,a_end,b_start,b_end,c_start,c_end)
f_inequalities(b_start,b_end,c_start,c_end,a_start,a_end)
f_inequalities(c_start,c_end,a_start,a_end,b_start,b_end)
            
def k_inequalities(a_start,a_end,b_start,b_end,c_start,c_end):
    for i in range(a_start + 1, a_end):  
        inequality = ""        
        inequality = inequality + f' + FK({i-1},{i})'
        
        for j in range(b_start + 1, b_end):  
            inequality = inequality + f' + KF({i},{j})'
            inequality = inequality + f' + KK({i},{j})'
            
        for j in range(c_start + 1, c_end):  
            inequality = inequality + f' + FK({j},{i})'
            inequality = inequality + f' + KK({j},{i})'
        
        inequality = inequality + ' <= 1\n'
        OUT.write (inequality)

k_inequalities(a_start,a_end,b_start,b_end,c_start,c_end)
k_inequalities(b_start,b_end,c_start,c_end,a_start,a_end)
k_inequalities(c_start,c_end,a_start,a_end,b_start,b_end)

# def kf_inequalities_1(a_start,a_end,b_start,b_end):
#     for i in range(a_start + 1, a_end):  
#         inequality = ""        
#         for j in range(b_start + 1, b_end):  
#             inequality = inequality + f' + KF({i},{j})'
#             inequality = inequality + f' + FK({i-1},{i})'        
#             inequality = inequality + ' <= 1\n'
#             OUT.write (inequality)
            
# def kf_inequalities_2(a_start,a_end,b_start,b_end):
#     for i in range(a_start + 1, a_end):  
#         inequality = ""        
#         for j in range(b_start + 1, b_end):  
#             inequality = inequality + f' + KF({i},{j})'
#             inequality = inequality + f' + FK({j},{j+1})'        
#             inequality = inequality + ' <= 1\n'
#             OUT.write (inequality)
            
# kf_inequalities_1(a_start,a_end,b_start,b_end)
# kf_inequalities_1(b_start,b_end,c_start,c_end)
# kf_inequalities_1(c_start,c_end,a_start,a_end)
# kf_inequalities_2(a_start,a_end,b_start,b_end)
# kf_inequalities_2(b_start,b_end,c_start,c_end)
# kf_inequalities_2(c_start,c_end,a_start,a_end)
# def fk_inequalities_1(a_start,a_end,b_start,b_end,c_start,c_end):
#     for i in range(a_start + 1, a_end):  
#             inequality = ""        
#             for j in range(b_start + 1, b_end):  
#                 inequality = inequality + f' + FK({i},{j})'
#                 inequality = inequality + f' + FK({i},{i+1})'
                
#             for j in range(c_start + 1, c_end):          
#                 # inequality = inequality + f' + FF({j},{i})'
#                 inequality = inequality + f' + FK({j},{i})'
#                 inequality = inequality + f' + FK({j},{j+1})'
#             inequality = inequality + ' <= 1'
#             OUT.write (inequality)
#             OUT.write ("\n\n")

# fk_inequalities_1(a_start,a_end,b_start,b_end,c_start,c_end)
# fk_inequalities_1(b_start,b_end,c_start,c_end,a_start,a_end)
# fk_inequalities_1(c_start,c_end,a_start,a_end,b_start,b_end)

# def fk_inequalities_2(a_start,a_end,b_start,b_end,c_start,c_end):
#     for i in range(a_start + 1, a_end):  
#             inequality = ""        
#             for j in range(b_start + 1, b_end):          
#                 # inequality = inequality + f' + FF({i},{j})'
#                 inequality = inequality + f' + FK({i},{j})'
#                 inequality = inequality + f' + FK({j},{j-1})'
                
#             for j in range(c_start + 1, c_end):          
#                 # inequality = inequality + f' + FF({j},{i})'
#                 inequality = inequality + f' + FK({j},{i})'
#                 inequality = inequality + f' + FK({i-1},{i})'
#             inequality = inequality + ' <= 1'
#             OUT.write (inequality)
#             OUT.write ("\n\n")

# fk_inequalities_2(a_start,a_end,b_start,b_end,c_start,c_end)
# fk_inequalities_2(b_start,b_end,c_start,c_end,a_start,a_end)
# fk_inequalities_2(c_start,c_end,a_start,a_end,b_start,b_end)

# stackings inequalities for max one stacking for k(5') face 
# def kf_inequalities_1(a_start,a_end,b_start,b_end,c_start,c_end):
#     for i in range(a_start + 1, a_end):  
#             inequality = ""        
#             for j in range(b_start + 1, b_end):
#                 inequality = inequality + f' + KF({i},{j})'
#                 # inequality = inequality + f' + KK({i},{j})'
#                 inequality = inequality + f' + FK({i-1},{i})'
                
#             for j in range(c_start + 1, c_end):
#                 inequality = inequality + f' + KF({j},{i})'
#                 # inequality = inequality + f' + KK({j},{i})'
#                 inequality = inequality + f' + FK({j-1},{j})'
#             inequality = inequality + ' <= 1'
#             OUT.write (inequality)
#             OUT.write ("\n\n")

# kf_inequalities_1(a_start,a_end,b_start,b_end,c_start,c_end)
# kf_inequalities_1(b_start,b_end,c_start,c_end,a_start,a_end)
# kf_inequalities_1(c_start,c_end,a_start,a_end,b_start,b_end)

# def kf_inequalities_2(a_start,a_end,b_start,b_end,c_start,c_end):
#     for i in range(a_start + 1, a_end):  
#             inequality = ""        
#             for j in range(b_start + 1, b_end):
#                 inequality = inequality + f' + KF({i},{j})'
#                 # inequality = inequality + f' + KK({i},{j})'
#                 inequality = inequality + f' + FK({j},{j+1})'
                
#             for j in range(c_start + 1, c_end):
#                 inequality = inequality + f' + KF({j},{i})'
#                 # inequality = inequality + f' + KK({j},{i})'
#                 inequality = inequality + f' + FK({i},{i+1})'
#             inequality = inequality + ' <= 1'
#             OUT.write (inequality)
#             OUT.write ("\n\n")

# kf_inequalities_2(a_start,a_end,b_start,b_end,c_start,c_end)
# kf_inequalities_2(b_start,b_end,c_start,c_end,a_start,a_end)
# kf_inequalities_2(c_start,c_end,a_start,a_end,b_start,b_end)
                
OUT.write("\n\n")

# OUT.write ("E <= " + str(loop8))
# OUT.write("\n")
OUT.write (listS + " <= " + str(loop8_s))
OUT.write("\n")
# OUT.write (listQ + " >= " + str(loop8_q))
# OUT.write("\n")
# OUT.write (listP + " <= " + str(loop8_p))
# OUT.write("\n")
OUT.write (listP + listQ + " <= " + str(loop8_pq))
OUT.write("\n")


def hard_constraints():
    
     OUT.write (f'WW({a_end},{b_start}) = 1')
     OUT.write ("\n\n")
     OUT.write (f'WW({b_end},{c_start}) = 1')
     OUT.write ("\n\n")
     OUT.write (f'WW({c_end},{a_start}) = 1')
     OUT.write ("\n\n")
     
hard_constraints()
            
OUT.write("\n")
OUT.write ("binary \n")

def binaries(a_start,a_end,b_start,b_end):
    for i in range(a_start + 1, a_end):
        for j in range(b_start + 1, b_end):           
            if RNA[i-1] + RNA[j-1] not in NP_list:
                # print(f'{i},{j}: {RNA[i-1]}{RNA[j-1]}')
                OUT.write (f'WW({i},{j}) \n')
            OUT.write (f'HH({i},{j}) \n')
            OUT.write (f'SS({i},{j}) \n')
            OUT.write (f'WH({i},{j}) \n')
            OUT.write (f'WH({j},{i}) \n')
            OUT.write (f'WS({i},{j}) \n')
            OUT.write (f'WS({j},{i}) \n')
            OUT.write (f'HS({i},{j}) \n')
            OUT.write (f'HS({j},{i}) \n')

binaries(a_start,a_end,b_start,b_end)
binaries(b_start,b_end,c_start,c_end)
binaries(c_start,c_end,a_start,a_end)

def q_binaries(a_start,a_end,b_start,b_end):
    for i in range(a_start+1,a_end):
        for j in range(b_start+1,b_end):
            OUT.write (f'FF({i},{j}) \n')
            OUT.write (f'FK({i},{j}) \n')
            OUT.write (f'KF({i},{j}) \n')
            OUT.write (f'KK({i},{j}) \n')
            
q_binaries(a_start,a_end,b_start,b_end)
q_binaries(b_start,b_end,c_start,c_end)
q_binaries(c_start,c_end,a_start,a_end)

def binaries_fk(a_start,a_end):
    for i in range(a_start, a_end):
        OUT.write (f'FK({i},{i+1}) \n')
        # print(f'FK({i},{i+1}) \n')
        
binaries_fk(a_start,a_end)
binaries_fk(b_start,b_end)
binaries_fk(c_start,c_end)
    
                    
OUT.write("\n")
OUT.write ("general\n")
OUT.write("E\n")
OUT.write  ("end \n")
        
OUT.close()






