# from constants_paths import *
# from prepro_utils import *

# # find files that do not corresponding .ct or .seq paired file
# directory_path = archive_path  # Replace with your directory path

# # print("Files without corresponding pairs:")
# # for file in unpaired_files:
# #     ic(file)

# # create 'noctfiles' folder and move .seq files without .ct to it 
# source_directory = archive_path  # Replace with the path to your source directory
# destination_directory = archive_path  # Replace with the path to your destination directory
# new_folder_name = 'noctfiles'  # Replace with your desired new folder name


# # create folders to write .lp model and optimization .sol results
# lp_dir = create_dir(lp_path, lp_folder_name)
# sol_dir = create_dir(lp_path, sol_folder_name)
# dot_bracket_dir = create_dir(lp_path, dot_bracket_folder_name)
# dot_bracket_archive_dir = create_dir(archive_path, dot_bracket_archive_folder_name)
# rnastructure_fold_dir =  create_dir(archive_path, rnastructure_folder_name)
# dot_bracket_rnastructure_dir = create_dir(archive_path, dot_bracket_rnastructure_folder_name)

# # # get all paths to .seq files 
# # seq_dir = os.path.join(archive_path, seq_len_dir)
# # seq_files = get_filenames(seq_dir, '.seq')

# # # get all paths to rnastructure .ct files 
# # ct_rnastruct_dir = os.path.join(archive_path, rnastructure_fold_dir)
# # ct_rnastruct_files = get_filenames(ct_rnastruct_dir, '.ct')

# # path to sequence from archive ii os.path.join(dest_dir, new_folder_name)
# chain_dir = os.path.join(archive_path, seq_len_dir)
# seq_files = get_filenames(chain_dir, '.seq')

# ic(seq_files[seq_number])
# ic(len(seq_files))

# chain_file = seq_files[seq_number]
# chain_name_with_ext = os.path.basename(chain_file)
# chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
# import os
# from pathlib import Path

# import subprocess

# from prepro_utils import *
# from prepro_run import *

# from constants_paths import *


 # get all paths to .seq files 
# seq_dir = os.path.join(archive_path, seq_len_dir)
# seq_files = get_filenames(seq_dir, '.seq')
# rnastructure_fold_dir =  create_dir(archive_path, rnastructure_folder_name)
# change working directory to path to RNAstructure exe files
# rnastruct_path = "C:/Program Files/RNAstructure6.5/exe"
# os.chdir(rnastruct_path)

# rnastruct_fold(seq_files, 51, 52, rnastructure_fold_dir)
# ct2dot(ct_rnastruct_files, 51, 52, dot_bracket_rnastructure_dir)
