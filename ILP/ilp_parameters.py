import os

# sequence path
chain_dir = 'data'
chain_f = '2ke6_a'

chain_fs = []

chain_path = os.path.join(chain_dir, chain_f)

with open( chain_path, 'r') as file:
    RNA = file.read().rstrip()    
    print (RNA)

# distance parameters
minD = 3 # min allowed distance between paired nts, also min hairpin loop size
maxH = 35 # max hairpin loop size
noH = 30
minI = 1 # min internal loop unpaired region size
maxI = 6 # max internal loop one side size
B = 1 # one nt bulge loop
maxB = 7 # max bulge loop

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