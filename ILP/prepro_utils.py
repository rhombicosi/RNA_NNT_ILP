import os
# from pathlib import Path
from collections import defaultdict
from icecream import ic
import shutil
import subprocess

def get_filenames(dir_path, f_type):
    f_list = []
    # list all .seq files in the archive directory
    for file in os.listdir(dir_path):
        # check the files which are end with specific extension
        if file.endswith(f_type):
            # add path name of a selected file to the list of files
            f_list.append(os.path.join(dir_path, file))
            # print(os.path.join(dir_path, file))
        
    return f_list

def parse_ct_file(filename):
    """
    Parses a .ct file and returns a list of nucleotides with their pairings.

    Args:
    filename (str): The path to the .ct file.

    Returns:
    list of dict: Each dictionary contains information about a nucleotide.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    # First line contains the number of nucleotides and possibly a comment
    num_nucleotides = int(lines[0].split()[0])
    
    # List to store nucleotide information
    nucleotides = []
    
    # Parse each nucleotide line
    for line in lines[1:num_nucleotides + 1]:
        parts = line.split()
        nucleotide_info = {
            'number': int(parts[0]),
            'nucleotide': parts[1],
            'prev': int(parts[2]),
            'next': int(parts[3]),
            'pair': int(parts[4])
        }
        nucleotides.append(nucleotide_info)
    
    return nucleotides

def parse_seq_file(filename):
    """
    Parses a .seq file and returns a dictionary with the sequence identifier and sequence.

    Args:
    filename (str): The path to the .seq file.

    Returns:
    dict: A dictionary containing the sequence identifier and sequence.
    """
    sequence_data = {}

    with open(filename, 'r') as file:
        lines = file.readlines()

    # for line in lines:
    #     line = line.strip()
    #     print(line)
    
    sequence_data['identifier'] = lines[1]
    seq = ''.join(lines[2])
    sequence_data['sequence'] = seq[:-2]
    return sequence_data

def find_files_without_pairs(directory):
    # Create a dictionary to store files by their base name
    files_dict = defaultdict(set)
    
    # List all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.seq') or filename.endswith('.ct'):
            # Extract the base name and extension
            base_name, ext = os.path.splitext(filename)
            # Add the extension to the set of this base name
            files_dict[base_name].add(ext)
    
    # Find files without pairs
    files_without_pairs = []
    for base_name, extensions in files_dict.items():
        if len(extensions) == 1:
            # There is only one type of file for this base name
            ext = list(extensions)[0]
            files_without_pairs.append(base_name + ext)
    
    return files_without_pairs

def create_directory_and_move_selected_files(src_dir, dest_dir, new_folder_name, selected_files):
    # Create the new directory
    new_folder_path = os.path.join(dest_dir, new_folder_name)
    os.makedirs(new_folder_path, exist_ok=True)
    
    # Move selected files from the source directory to the new directory
    for filename in selected_files:
        src_file = os.path.join(src_dir, filename)
        
        # Check if the file exists in the source directory
        if os.path.isfile(src_file):
            # Move the file to the new directory
            shutil.move(src_file, new_folder_path)
            # ic(f"Moved: {filename}")
        else:
            ic(f"File not found: {filename}")
    
    ic(f"Selected files have been moved to {new_folder_path}")

def create_directory_and_copy_selected_files(src_dir, dest_dir, new_folder_name, selected_files):
    # Create the new directory
    new_folder_path = os.path.join(dest_dir, new_folder_name)
    os.makedirs(new_folder_path, exist_ok=True)
    
    # Copy selected files from the source directory to the new directory
    for filename in selected_files:
        src_file = os.path.join(src_dir, filename)
        
        # Check if the file exists in the source directory
        if os.path.isfile(src_file):
            # Copy the file to the new directory
            shutil.copy(src_file, new_folder_path)
            # ic(f"Copied: {filename}")
        else:
            print(f"File not found: {filename}")
    
    print(f"Selected files have been copied to {new_folder_path}")

# select sequences of a specific length less then some value
def get_seq_of_len(seq_list, seq_len, ct_list):
    seq_len_list = []
    ct_len_list = []
    seq_len_files = []
    ct_len_files = []

    for (seq,ct) in zip(seq_list, ct_list):
        seq_data = parse_seq_file(seq)
        ct_data = parse_ct_file(ct)        

        if len(seq_data['sequence']) <= seq_len:
            seq_len_files.append(seq)
            ct_len_files.append(ct)
            seq_len_list.append(seq_data)
            ct_len_list.append(ct_data)
            # ic(seq)
            # ic(f"Identifier: {seq_data['identifier']}")
            # ic(f"Sequence: {seq_data['sequence']}")

    # ic(seq_len_list)
    # ic(ct_len_list)

    return seq_len_files,ct_len_files

# new_dir_name = 'lp'; rel_path_to_save = '../'
def create_dir(path_to_save, new_dir_name):
    # file_parent_dir = Path(__file__).parent
    # path_to_save = (file_parent_dir/rel_path_to_save).resolve()
    new_dir_path = os.path.join(path_to_save, new_dir_name)
    os.makedirs(new_dir_path, exist_ok=True)

    return new_dir_path

# convert .ct files to dot-bracket notation
def ct2dot(ct_files, first_file, last_file, dot_bracket_dir):
    for ct_number in range(first_file, last_file):
        ct = ct_files[ct_number]
        ct_name_with_ext = os.path.basename(ct)
        ct_name_without_ext = os.path.splitext(ct_name_with_ext)[0]

        result = subprocess.run(['ct2dot', ct,'1', f'{dot_bracket_dir}/{ct_name_without_ext + ".txt"}'], capture_output=True, text=True)
        print(result.stdout)

# generate .ct with RNA structure
def rnastruct_fold(seq_files, first_file, last_file, fold_dir):
    for seq_number in range(first_file, last_file):
        seq = seq_files[seq_number]
        seq_name_with_ext = os.path.basename(seq)
        seq_name_without_ext = os.path.splitext(seq_name_with_ext)[0]

        result = subprocess.run(['fold', seq, f'{fold_dir}/{seq_name_without_ext + ".ct"}'], capture_output=True, text=True)
        print(result.stdout)


def dot_from_txt(f_txt):
    # Open the file in read mode 'your_file.txt'
    with open(f_txt, 'r') as file:
        # Read all lines into a list
        lines = file.readlines()
        
        # Check if the file has at least three lines
        if len(lines) >= 3:
            # Access the third line (index 2 because indexing starts from 0)
            third_line = str(lines[2]).strip()
            return(third_line)
        else:
            pass
            # print("The file does not have three lines.")
