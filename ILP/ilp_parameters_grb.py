import os
from icecream import ic

# sequence path
chain_dir = 'data'
chain_f = '2ku0_a' #'2k5z'#'2KE6_A' 

chain_fs = []

chain_path = os.path.join(chain_dir, chain_f)

with open( chain_path, 'r') as file:
    RNA = file.read().rstrip() 
    sizeRNA = len(RNA)
    ic(chain_f,RNA)
    ic(sizeRNA)


# MFE parameter //  reference structure energy by RNAEval
MFE = -970


# distance parameters
minD = 3 # min allowed distance between paired nts, also min hairpin loop size
maxH = 35 # max hairpin loop size
noH = 30 # initiation params for hairpins bigger than noH are approximated
minI = 1 # min internal loop unpaired region size
maxI = 35 # max internal loop one side size
noI = 30 # initiation params for internals bigger than noH are approximated
maxB = 35 # max bulge loop
noB = 30

# number of loops parameters
numH = len(RNA)//5
numI = len(RNA)//4
numB = len(RNA)//3

# print(numH)
# print(numI)
# print(numB)

# R = 1.9872036 × 10-3	kcal.K-1.mol-1 is the gas constant and T is the absolute temperature, 310.15 K
RT = 0.616
Cbulge = -0.9