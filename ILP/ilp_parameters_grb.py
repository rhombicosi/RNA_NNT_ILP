# fixed params for creating ILP model

# MFE parameter //  reference structure energy by RNAEval
MFE = -1500

# distance parameters
minD = 3 # min allowed distance between paired nts, also min hairpin loop size
maxH = 9 # max hairpin loop size
noH = 9 # initiation params for hairpins bigger than noH are approximated
minI = 1 # min internal loop unpaired region size
maxI = 6 # max internal loop one side size
noI = 6 # initiation params for internals bigger than noH are approximated #TODO 22 is max size of one side so max loop is 44
maxB = 5 # max bulge loop
noB = 5
maxM = 5 # max size of one side of the multiloop
noM = 15 # total max size of the internal loop
L = 15 # maximum number of unpaired nucleotides

# R = 1.9872036 × 10-3	kcal.K-1.mol-1 is the gas constant and T is the absolute temperature, 310.15 K
RT = 0.616
Cbulge = -0.9

# multiloop params
a = 3.4 # penalty for closing the multiloop
b = 0.9 # penalty per branch (base pair)
c = 0.4 # penalty per unpaired nucleotide