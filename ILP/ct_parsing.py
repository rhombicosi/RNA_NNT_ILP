import os
from icecream import ic

from pathlib import Path

# sequence path
# chain_dir = 'data'
# chain_f = '2ku0_a' #'2k5z'#'2KE6_A' 

# chain_fs = []
cwd = Path.cwd()
mod_path = Path(__file__).parent.parent

relative_path = '../../ARCHIVE II/'
src_path = (mod_path /relative_path).resolve()


filename_ct = os.path.join(src_path, '5s_Acanthamoeba-castellanii-1.ct')


print(filename_ct)

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

# Example usage
# filename = 'example.ct'
nucleotide_data = parse_ct_file(filename_ct)

# Print parsed data
for nucleotide in nucleotide_data:
    print(nucleotide)

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

# Example usage
filename_seq = os.path.join(src_path, '5s_Acanthamoeba-castellanii-1.seq')
seq_data = parse_seq_file(filename_seq)

# Print parsed data
print(f"Identifier: {seq_data['identifier']}")
print(f"Sequence: {seq_data['sequence']}")
