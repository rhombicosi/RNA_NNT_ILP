from pathlib import Path
import pandas as pd

# archive ii path
# with .ct and .seq files
cwd = Path.cwd()
code_path = Path(__file__).parent.parent
arch_rel_path = '../../ARCHIVE II/'
archive_path = (code_path/arch_rel_path).resolve()

# folder name for files that should be excluded from test
badfiles_dir = 'noctfiles'  

# length for sequences to be tested
seq_len = 60

# folders to save .seq and .ct files with sequences of seq_len
seq_len_dir = f'RNA_seq_{seq_len}'
ct_len_dir = f'RNA_ct_{seq_len}'

# folders to save lp models and solutions
ilp_parent_dir = Path(__file__).parent
lp_rel_path = '../'
lp_path = (ilp_parent_dir/lp_rel_path).resolve()
lp_folder_name = 'lp_multi'
sol_folder_name = 'sol_multi'
solstart_folder_name = 'sol_start'
dot_bracket_folder_name = f'dot_bracket_{seq_len}'
dot_bracket_archive_folder_name = f'dot_bracket_archive_{seq_len}'
rnastructure_folder_name = f'RNAstructure_fold_{seq_len}'
dot_bracket_rnastructure_folder_name = f'dot_bracket_RNAstructure_{seq_len}'
efn2_archive_folder_name = f'efn2_archive_{seq_len}'

# RNAstructure path
rnastruct_path = "C:/Program Files/RNAstructure6.5/exe"

# index of the sequence from archive ii currently under test
# seq_number = 51

# objective val
obj_val = 0

# store generated .lp files
lp_files = []
