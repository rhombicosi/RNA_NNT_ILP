import re
from binary_variables_grb import *
from constants_paths import *
from prepro_run import *
from prepro_utils import *

def pairs2brackets(filepath, RNA): 
    lngth = len(RNA)
    print(lngth)

    pattern = "\((.*?)\)"

    fold = ["." for _ in range(lngth)] 
        
    with open(filepath) as fp:

        line = fp.readline()      
        cnt = 1
        pairs = []
        
        while line: 
            
            if " 1" in line and "Q(" in line:
                print("{}".format(line.strip()))
            
            if " 1" in line and "F(" in line:
                print("{}".format(line.strip()))

            if " 1" in line and "L(" in line:
                print("{}".format(line.strip()))

            if " 1" in line and "H(" in line:
                print("{}".format(line.strip()))
            
            if " 1" in line and "I(" in line:
                print("{}".format(line.strip()))

            if " 1" in line and "B(" in line:
                print("{}".format(line.strip()))

            if " 1" in line and "P(" in line:
                print("{}".format(line.strip()))
                substring = re.search(pattern, line).group(1)
                indices = substring.split(",")
            
                i = int(indices[0])
                j = int(indices[1])

                pairs.append((i,j))
                fold[i - 1] = "("
                fold[j - 1] = ")"
                
            line = fp.readline()
            # p_str = line.find("P(")
            cnt += 1

    fold = ''.join([str(elem) for elem in fold])

    return(fold,pairs)

def brackets2pairs(dot_bracket):
    pair_stack = []
    base_pairs = []

    for i, symbol in enumerate(dot_bracket, start=1):
        if symbol in "({[":
            pair_stack.append(i)
        elif symbol in ")}]":
            if pair_stack:
                j = pair_stack.pop()
                base_pairs.append((j, i))

    sorted_bp = sorted(base_pairs, key=lambda x: x[0])

    return sorted_bp

def compare2folds(generated, reference):
    beta = 2
    TP = len(reference.intersection(generated)) # pairs that have been identified correctly
    FN = len(reference - generated) # pairs that have not been identified
    FP = len(generated - reference) # pairs that have been identified as part of the structure incorrectly

    # precision = TP/(TP + FP)
    # recall = TP/(TP + FN)
    f1 = 2*TP/(2*TP + FP + FN)
    fbeta = (1+beta**2)*TP/((1+beta**2)*TP + (beta**2)*FN + FP)

    # print(TP)
    # print(FN)
    # print(FP)

    return (f1,fbeta)

def sol_analyse(seq_files, seq_number, sol_dir, dot_bracket_archive_dir, dot_bracket_rnastructure_dir):

    chain_file = seq_files[seq_number]
    chain_name_with_ext = os.path.basename(chain_file)        
    chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
    lp_file_name = chain_name_without_ext

    seq_data = parse_seq_file(chain_file)
    this_RNA = seq_data['sequence']

    filepath = os.path.join(sol_dir, f'{lp_file_name}-loopdeco.sol')

    ic(filepath)

    # archive referece
    ref_bracket_path = f'{dot_bracket_archive_dir}/{lp_file_name + ".txt"}'
    ref_brackets = dot_from_txt(ref_bracket_path)

    # rnastructure reference
    rnastruct_bracket_path = f'{dot_bracket_rnastructure_dir}/{lp_file_name + ".txt"}'
    rnastruct_brackets = dot_from_txt(rnastruct_bracket_path)


    (gen_brackets, gen_pairs) = pairs2brackets(filepath, this_RNA)

    print(gen_brackets)
    print(gen_pairs)

    # pairs generated by the model and pairs of a reference structure
    generated = set(gen_pairs)
    reference = set(brackets2pairs(ref_brackets))
    rnastruct = set(brackets2pairs(rnastruct_brackets))

    print(generated)
    print(reference)
    print(rnastruct)

    (f1_ref,fbeta_ref) = compare2folds(generated, reference)
    (f1_rnastruct,fbeta_rnastruct) = compare2folds(generated, rnastruct)


    print(f1_ref)
    print(fbeta_ref)

    print(f1_rnastruct)
    print(fbeta_rnastruct)

    return f1_ref, fbeta_ref, f1_rnastruct, fbeta_rnastruct