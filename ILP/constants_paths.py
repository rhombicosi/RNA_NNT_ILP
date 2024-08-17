from pathlib import Path

# archive ii path
# with .ct and .seq files
cwd = Path.cwd()
code_path = Path(__file__).parent.parent
arch_rel_path = '../../ARCHIVE II/'
archive_path = (code_path/arch_rel_path).resolve()

# length for sequences to be tested
seq_len = 60

# folders to save .seq and .ct files with sequences of seq_len
seq_len_dir = f'RNA_seq_{seq_len}'
ct_len_dir = f'RNA_ct_{seq_len}'

# folders to save lp models and solutions
ilp_parent_dir = Path(__file__).parent
lp_rel_path = '../'
lp_path = (ilp_parent_dir/lp_rel_path).resolve()
lp_folder_name = 'lp'
sol_folder_name = 'sol'
dot_bracket_folder_name = f'dot_bracket_{seq_len}'
dot_bracket_archive_folder_name = f'dot_bracket_archive_{seq_len}'
rnastructure_folder_name = f'RNAstructure_fold_{seq_len}'
dot_bracket_rnastructure_folder_name = f'dot_bracket_RNAstructure_{seq_len}'

# index of the sequence from archive ii currently under test
seq_number = 1
