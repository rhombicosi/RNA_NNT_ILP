import os
from icecream import ic
from archive_prepro_parse import *

# path to archive ii .ct and .seq files 
# t = archive_path
# seq_len_dir = f'RNA_seq_{seq_len}'
# ct_len_dir = f'RNA_ct_{seq_len}'
# sequences = get_filenames(archive_path, ".seq")
# ct_list = get_filenames(archive_path, ".ct")



# sequence path os.path.join(dest_dir, new_folder_name)
chain_dir = os.path.join(archive_path, seq_len_dir)
seq_files = get_filenames(chain_dir, '.seq')

seq_number = 1
ic(seq_files[seq_number])
ic(len(seq_files))

chain_path = seq_files[seq_number]
chain_f = f'seq-{seq_number}' #os.path.basename(chain_path)

# with open(chain_path, 'r') as file:
seq_data = parse_seq_file(chain_path)
RNA = seq_data['sequence']
# RNA = file.read().rstrip() 
sizeRNA = len(RNA)
ic(chain_f,RNA)
ic(sizeRNA)



with open( chain_path, 'r') as file:

    RNA = seq_data['sequence']
    sizeRNA = len(RNA)
    ic(chain_f,RNA)
    ic(sizeRNA)


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