from fnmatch import fnmatchcase
import math
            
from score_matrix import bp_df, stack_df 

sqnc_name = "1etg"

with open(sqnc_name, 'r') as file:
    MOTIF = file.read().rstrip()    
    # print (MOTIF)

length = len(MOTIF)
# print (length)

SCORE = 116880

length_for_range = length + 1
weights = [[0 for x in range(1, length_for_range + 1)] for y in range(1, length_for_range + 1)]

# strands
a_start = 1
a_end = 4

b_start = a_end + 1
b_end = length

OUT = open(f'{sqnc_name}-motifs-dl.lp', "w")

OUT.write ("Maximize \n")
OUT.write ("E \n")

listP = "" # listP will sum the edge-to-edge base pairs variables

# list of canonical pairs
CP_list = ['AU','UA','CG','GC','GU','UG']

# all edge2edge non canonical contacts
def listp_construct(a_start,a_end,b_start,b_end,listP):
    
    for i in range(a_start, a_end+1):
        for j in range(b_start, b_end+1):           
                
            # edge2edge contacts weights
            ww_weight = bp_df.loc['WW_cis', MOTIF[i-1] + MOTIF[j-1]]
            hh_weight = bp_df.loc['HH_cis', MOTIF[i-1] + MOTIF[j-1]]
            ss_weight = bp_df.loc['SS_cis', MOTIF[i-1] + MOTIF[j-1]]
            
            ws_weight = bp_df.loc['WS_cis', MOTIF[i-1] + MOTIF[j-1]]
            wh_weight = bp_df.loc['WH_tran', MOTIF[i-1] + MOTIF[j-1]]
            hs_weight = bp_df.loc['HS_tran', MOTIF[i-1] + MOTIF[j-1]]
            
            if not math.isnan(float(ww_weight)):                
                if MOTIF[i-1] + MOTIF[j-1] not in CP_list:            
                    listP = listP + f' + {ww_weight*0.5} WW({i},{j})'
                    listP = listP + f' + {ww_weight*0.5} WW({j},{i})'
                if (i == a_start and j == b_end) or (i == a_end and j == b_start):
                    listP = listP + f' + {ww_weight*0.5} WW({i},{j})'
                    listP = listP + f' + {ww_weight*0.5} WW({j},{i})'
                
            if not math.isnan(float(hh_weight)):
                listP = listP + f' + {hh_weight*0.5} HH({i},{j})'
                listP = listP + f' + {hh_weight*0.5} HH({j},{i})'
            if not math.isnan(float(ss_weight)):
                listP = listP + f' + {ss_weight*0.5} SS({i},{j})'
                listP = listP + f' + {ss_weight*0.5} SS({j},{i})'
                
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

# print(listP)

OUT.write ("\nsuch that \n")
OUT.write (listP + " - E = 0 \n\n")

# The following lines generate the inequalities to ensure
# that each edge is paired to at most one other edge

# WW edge inequalities
def ww_inequalities(a_start,a_end,b_start,b_end):
    for i in range(a_start, a_end+1):
        inequality = ""    
        for j in range(b_start, b_end+1):
            if MOTIF[i-1] + MOTIF[j-1] not in CP_list:
                inequality = inequality + f' + WW({i},{j})'
            if (i == a_start and j == b_end) or (i == a_end and j == b_start):
                inequality = inequality + f' + WW({i},{j})'

            inequality = inequality + f' + WH({i},{j})'
            inequality = inequality + f' + WS({i},{j})'
        
        inequality = inequality + ' <= 1'
        OUT.write (inequality)           		
        OUT.write ("\n\n")

ww_inequalities(a_start,a_end,b_start,b_end)
ww_inequalities(b_start,b_end,a_start,a_end)
        
# HH edge inequalities
def  hh_inequalities(a_start,a_end,b_start,b_end):
    for i in range(a_start, a_end+1):
        inequality = ""    
        for j in range(b_start, b_end+1):

            inequality = inequality + f' + HH({i},{j})'
            inequality = inequality + f' + WH({j},{i})'
            inequality = inequality + f' + HS({i},{j})'
                        
        inequality = inequality + ' <= 1'
        OUT.write (inequality)           		
        OUT.write ("\n\n")
            
hh_inequalities(a_start,a_end,b_start,b_end)
hh_inequalities(b_start,b_end,a_start,a_end)

# SS edge inequalities
def ss_inequalities(a_start,a_end,b_start,b_end):
    for i in range(a_start, a_end+1):
        inequality = ""    
        for j in range(b_start, b_end+1):
            
            inequality = inequality + f' + SS({i},{j})'        
            inequality = inequality + f' + WS({j},{i})'
            inequality = inequality + f' + HS({j},{i})'
                        
        inequality = inequality + ' <= 1'
        OUT.write (inequality)           		
        OUT.write ("\n\n")
            
        # print(f'inequality #{i}: {inequality}')
            
ss_inequalities(a_start,a_end,b_start,b_end)
ss_inequalities(b_start,b_end,a_start,a_end)

# equalities

def equalities(a_start,a_end,b_start,b_end):
    equality = ""    
    for i in range(a_start, a_end+1):
        for j in range(b_start, b_end+1):
            
            equality = f'WW({i},{j}) - WW({j},{i}) = 0' 
            OUT.write (equality)           		
            OUT.write ("\n\n")
            equality = f'HH({i},{j}) - HH({j},{i}) = 0' 
            OUT.write (equality)           		
            OUT.write ("\n\n")
            equality = f'SS({i},{j}) - SS({j},{i}) = 0' 
            OUT.write (equality)           		
            OUT.write ("\n\n") 

equalities(a_start,a_end,b_start,b_end)

def hard_constraints():
    
     OUT.write (f'WW({a_end},{b_start}) = 1')
     OUT.write ("\n\n")
     OUT.write (f'WW({b_end},{a_start}) = 1')
     OUT.write ("\n\n")
     
hard_constraints()

def none_score():    

    for i in range(a_start, a_end+1):
        for j in range(b_start, b_end+1): 
            ss_weight = bp_df.loc['SS_cis', MOTIF[i-1] + MOTIF[j-1]]
            if math.isnan(float(ss_weight)):
                OUT.write (f'SS({i},{j}) = 0\n')
                OUT.write (f'SS({j},{i}) = 0\n')

none_score()

OUT.write(f'E <= {SCORE} \n')
            
OUT.write("\n")
OUT.write ("binary \n")

def binaries(a_start,a_end,b_start,b_end):
    for i in range(a_start, a_end+1):
        for j in range(b_start, b_end+1):           
            if MOTIF[i-1] + MOTIF[j-1] not in CP_list:
                # print(f'{i},{j}: {RNA[i-1]}{RNA[j-1]}')
                OUT.write (f'WW({i},{j}) \n')
                OUT.write (f'WW({j},{i}) \n')
            if (i == a_start and j == b_end) or (i == a_end and j == b_start):
                OUT.write (f'WW({i},{j}) \n')
                OUT.write (f'WW({j},{i}) \n')
            OUT.write (f'HH({i},{j}) \n')
            OUT.write (f'HH({j},{i}) \n')
            OUT.write (f'SS({i},{j}) \n')
            OUT.write (f'SS({j},{i}) \n')
            OUT.write (f'WH({i},{j}) \n')
            OUT.write (f'WH({j},{i}) \n')
            OUT.write (f'WS({i},{j}) \n')
            OUT.write (f'WS({j},{i}) \n')
            OUT.write (f'HS({i},{j}) \n')
            OUT.write (f'HS({j},{i}) \n')

binaries(a_start,a_end,b_start,b_end)    
                    
OUT.write("\n")
OUT.write ("general\n")
OUT.write("E\n")
OUT.write  ("end \n")
        
OUT.close()






