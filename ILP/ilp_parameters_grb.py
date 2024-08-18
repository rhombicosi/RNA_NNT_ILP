from icecream import ic
from prepro_utils import *
from constants_paths import *
import prepro_run

# chain_f = f'seq-{prepro_run.seq_number}' 

# with open(chain_path, 'r') as file:
seq_data = parse_seq_file(prepro_run.chain_file)
RNA = seq_data['sequence']
# RNA = file.read().rstrip() 
sizeRNA = len(RNA)
# ic(chain_f,RNA)
ic(sizeRNA)

with open( prepro_run.chain_file, 'r') as file:

    RNA = seq_data['sequence']
    sizeRNA = len(RNA)
    # ic(chain_f,RNA)
    # ic(sizeRNA)


# MFE parameter //  reference structure energy by RNAEval
MFE = -970

# distance parameters
minD = 3 # min allowed distance between paired nts, also min hairpin loop size
maxH = 20 # max hairpin loop size
noH = 22 # initiation params for hairpins bigger than noH are approximated
minI = 1 # min internal loop unpaired region size
maxI = 22 # max internal loop one side size
noI = 22 # initiation params for internals bigger than noH are approximated
maxB = 20 # max bulge loop
noB = 20

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