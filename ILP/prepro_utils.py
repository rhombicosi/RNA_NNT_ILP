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
    
    # ic(f"Selected files have been moved to {new_folder_path}")

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
    
    # print(f"Selected files have been copied to {new_folder_path}")

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

        ct = os.path.join(fold_dir,f'{seq_name_without_ext + ".ct"}')

        result = subprocess.run(['fold', seq, ct], capture_output=True, text=True)
        print(result.stdout)

# calculate energy based on .ct file data
def rnastruct_efn2(ct_files, first_file, last_file, efn2_dir):
    for ct_number in range(first_file, last_file):
        ct = ct_files[ct_number]
        ct_name_with_ext = os.path.basename(ct)
        ct_name_without_ext = os.path.splitext(ct_name_with_ext)[0]

        efn2 = os.path.join(efn2_dir,f'{ct_name_without_ext + ".txt"}')

        result = subprocess.run(['efn2', ct, efn2], capture_output=True, text=True)
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

def get_energy_from_ct_file(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
        # Check if "ENERGY" is in the first line and extract the value
        if "ENERGY" or "Energy" in first_line:
            # Split the line by spaces and get the ENERGY value
            parts = first_line.split()
            for i, part in enumerate(parts):
                if part == "ENERGY" or part == "Energy":
                    # The value should be the next part after "ENERGY ="
                    energy_value = float(parts[i + 2])
                    return energy_value
        else:
            print("ENERGY value not found in the first line.")
            return None
        
def save_data_to_txt(file_name, data):
    # Define the headers
    headers = ["RNA sequence name", "# of nts", "MFE ILP", "MFE ARCHIVE", "MFE RNAstr",
               "F1 ILP", "F1 RNAstr", "Fb ILP", "Fb RNAstr"]

    # Check if the file exists and is not empty
    file_exists = os.path.exists(file_name) and os.path.getsize(file_name) > 0

    # Open the file for appending
    with open(file_name, 'a') as file:
        # If the file does not exist or is empty, write the headers
        if not file_exists:
            file.write("\t".join(headers) + "\n")
        
        # Write the data rows
        for row in data:
            # Convert all items in row to string and join with tab
            file.write("\t".join(map(str, row)) + "\n")

def write_results_to_file(sequence_name, rna_len, time, mfe_gen, mfe_ref, mfe_rna, f1_gen, f1_rna, fb_gen, fb_rna, mcc_gen, mcc_rna, filename="ilp_results.txt"):
    # Define the headers
    headers = ["RNA sequence name", "# of nts", "Time(s)", "MFE ILP", "MFE ARCHIVE", "MFE RNAstr",
               "F1 ILP", "F1 RNAstr", "Fb ILP", "Fb RNAstr", "INF ILP", "INF RNAstr"]
    

    # Check if the file exists and is not empty
    file_exists = os.path.exists(filename) and os.path.getsize(filename) > 0
    # Format the floating point numbers to 2 decimal places and ensure all values are strings
    values = [
        sequence_name, 
        rna_len,
        f"{time:.2f}",
        f"{mfe_gen:.2f}", 
        f"{mfe_ref:.2f}", 
        f"{mfe_rna:.2f}", 
        f"{f1_gen:.2f}", 
        f"{f1_rna:.2f}", 
        f"{fb_gen:.2f}", 
        f"{fb_rna:.2f}",
        f"{mcc_gen:.2f}",
        f"{mcc_rna:.2f}"

    ]

    headers_line = "{:<45}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\n".format(*headers)
    
    # Format the output so that each value is aligned with tabs
    line = "{:<45}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\n".format(*values)
    
    # Open the file in append mode, create if not exists
    with open(filename, 'a') as file:
        # If the file does not exist or is empty, write the headers
        if not file_exists:
            file.write(headers_line)
        
        # Write the data rows
        file.write(line)

def read_sol(f_name):
    solvars = {}
    with open(f_name,'r') as f:
            lines = f.readlines()
            for line in lines[2:]:
                parts = line.strip().split()
                if len(parts) >= 2:  # Ensure at least two parts exist
                    key, value = parts[0], parts[1]
                    solvars[key] = value
    
    return solvars
