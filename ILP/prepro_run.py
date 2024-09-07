import os
from prepro_utils import *
from constants_paths import *

# create folders to write .lp model and optimization .sol results
lp_dir = create_dir(lp_path, lp_folder_name)
sol_dir = create_dir(lp_path, sol_folder_name)
dot_bracket_dir = create_dir(lp_path, dot_bracket_folder_name)
dot_bracket_archive_dir = create_dir(archive_path, dot_bracket_archive_folder_name)
rnastructure_fold_dir =  create_dir(archive_path, rnastructure_folder_name)
dot_bracket_rnastructure_dir = create_dir(archive_path, dot_bracket_rnastructure_folder_name)
efn2_archive_dir = create_dir(archive_path, efn2_archive_folder_name)

# get .seq and .ct files after clearing redundant files
seq_list = get_filenames(archive_path, ".seq")
ct_list = get_filenames(archive_path, ".ct")

# create path to archive .ct files 
ct_archive_dir = os.path.join(archive_path, ct_len_dir)

# create path to rnastructure .ct files 
ct_rnastruct_dir = os.path.join(archive_path, rnastructure_fold_dir)

# create path to to .seq files 
seq_dir = os.path.join(archive_path, seq_len_dir)

if __name__ == "__main__":

    # find files that do not corresponding .ct or .seq paired file
    directory_path = archive_path  # Replace with your directory path
    unpaired_files = find_files_without_pairs(directory_path)

    # print("Files without corresponding pairs:")
    # for file in unpaired_files:
    #     ic(file)

    # create 'noctfiles' folder and move .seq files without .ct to it 
    source_directory = archive_path  # Replace with the path to your source directory
    destination_directory = archive_path  # Replace with the path to your destination directory
    new_folder_name = 'noctfiles'  # Replace with your desired new folder name
    selected_files = unpaired_files  # Replace with the list of files you want to move

    create_directory_and_move_selected_files(source_directory, destination_directory, new_folder_name, selected_files)

    # for z in zip(seq_list[0:1300],ct_list[0:1300]):
    #     ic(z)
    # ic(len(seq_list))
    # ic(len(ct_list))

    seq_len_files,ct_len_files = get_seq_of_len(seq_list,seq_len,ct_list)

    create_directory_and_copy_selected_files(source_directory, destination_directory, seq_len_dir, seq_len_files)
    create_directory_and_copy_selected_files(source_directory, destination_directory, ct_len_dir, ct_len_files)    

    # change working directory to path to RNAstructure exe files
    rnastruct_path = "C:/Program Files/RNAstructure6.5/exe"
    os.chdir(rnastruct_path)

    # get all paths to archive .ct files 
    ct_archive_files = get_filenames(ct_archive_dir, '.ct')

    # obtain dot-bracket notations for all archive .ct references
    ct2dot(ct_archive_files, 0, len(ct_archive_files), dot_bracket_archive_dir)

    # obtain enrgies for all archive .ct files 
    rnastruct_efn2(ct_archive_files, 0, len(ct_archive_files), efn2_archive_dir)

    # get all paths to .seq files    
    seq_files = get_filenames(seq_dir, '.seq')

    # generate secondary structures with RNAstructure algorithm, base settings
    rnastruct_fold(seq_files, 0, len(seq_files), rnastructure_fold_dir)   

    # get all paths to rnastructure .ct files
    ct_rnastruct_files = get_filenames(ct_rnastruct_dir, '.ct') 

    # obtain dot-bracket notations for all rnastructure .ct references
    ct2dot(ct_rnastruct_files, 0, len(ct_rnastruct_files), dot_bracket_rnastructure_dir)

# path to sequence from archive ii os.path.join(dest_dir, new_folder_name)
chain_dir = os.path.join(archive_path, seq_len_dir)
seq_files = get_filenames(chain_dir, '.seq')

ic(seq_files[seq_number])
ic(len(seq_files))

chain_file = seq_files[seq_number]
chain_name_with_ext = os.path.basename(chain_file)
chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]

print(f' # of sequences of lengths <= {seq_len} :: {len(seq_files)}')



# get all paths to efn2 .txt files 
# efn2_archive_files = get_filenames(efn2_archive_dir, '.txt')

# MFE_ref = get_energy_from_ct_file(efn2_archive_files[seq_number])
# MFE_rna = get_energy_from_ct_file(ct_rnastruct_files[seq_number])

# print(MFE_ref)
# print(MFE_rna)



