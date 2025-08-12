import os
from pathlib import Path
from lilp import *

seq_len = 60
seq_no = 46

cwd = Path.cwd()
code_path = Path(__file__).parent.parent
arch_rel_path = '../../ARCHIVE II/'
archive_path = (code_path/arch_rel_path).resolve()
seq_len_dir = f'RNA_seq_{seq_len}'

chain_dir = os.path.join(archive_path, seq_len_dir)
seq_files = get_filenames(chain_dir, '.seq')

chain_file = seq_files[seq_no]
chain_name_with_ext = os.path.basename(chain_file)        
chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
lp_file_name = chain_name_without_ext
seq_data = parse_seq_file(chain_file)

rna = seq_data['sequence']
print(rna)

bp1 = BasePair(3,54,rna)
bp2 = BasePair(12,53,rna)

# rna_model = LILP(rna)
# rna_model.create_variables(1, 1, 1, 1, 1)

# counter = 0


# for b in rna_model.bulge_loops:   
#     if b.var.VarName == 'BULGE_3_54_12_53':
#         print(b.size)
#     # print(f"[{counter}]:: Variable name: {b.var.VarName}, Energy: {b.energy}, Size: {b.size} Gurobi var: {b.var}")
#     counter+=1

# script_dir = os.path.dirname(os.path.abspath(__file__))
# parent_dir = os.path.dirname(script_dir)
# filepath = os.path.join(parent_dir, f'lilp-start-{seq_no}.sol')

# pairs2brackets(filepath, rna)
# calculate_sol_energy(filepath, rna)
# for f in range(15):
#     # filepath = os.path.join(parent_dir, f'lilp-start-{seq_no}-incumbent_{f}.sol')
#     filepath = os.path.join(parent_dir, f'lilp-start-{seq_no}.sol')
#     pairs2brackets(filepath, rna)