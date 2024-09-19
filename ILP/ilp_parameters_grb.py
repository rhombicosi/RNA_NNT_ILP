# fixed params for creating ILP model

# MFE parameter //  reference structure energy by RNAEval
MFE = -1500

# distance parameters
minD = 3 # min allowed distance between paired nts, also min hairpin loop size
maxH = 27 # max hairpin loop size
noH = 27 # initiation params for hairpins bigger than noH are approximated
minI = 1 # min internal loop unpaired region size
maxI = 22 # max internal loop one side size
noI = 22 # initiation params for internals bigger than noH are approximated
maxB = 20 # max bulge loop
noB = 20

# R = 1.9872036 × 10-3	kcal.K-1.mol-1 is the gas constant and T is the absolute temperature, 310.15 K
RT = 0.616
Cbulge = -0.9

# multiloop params
a = 3.4 # penalty for closing the multiloop
b = 0.9 # penalty per branch (base pair)
c = 0.4 # penalty per unpaired nucleotide