import pandas as pd
from utils.sol_converter import *
from utils.prepro_run import *
from utils.constants_paths import *
from optimize import *

def add_column(df, column_name, values):
    df[column_name] = values
    return df

ref_mfe_files = get_filenames(efn2_archive_dir, '.txt')
ref_mfe_names = [os.path.splitext(os.path.basename(f))[0] for f in ref_mfe_files]
results_df = pd.DataFrame(index=ref_mfe_names)

ref_MFEs = [get_energy_from_ct_file(f) for f in ref_mfe_files]

rna_mfe_names = get_filenames(ct_rnastruct_dir, '.ct')
rna_MFEs = [get_energy_from_ct_file(f) for f in rna_mfe_names]

vienna_mfe_names = get_filenames(dot_bracket_viennaRNA_dir, '_enrg.txt')
vienna_MFEs = [float(open(f).read().strip()) for f in vienna_mfe_names] 

add_column(results_df, 'MFE_ref', ref_MFEs)
add_column(results_df, 'MFE_rna', rna_MFEs)
add_column(results_df, 'MFE_vienna', vienna_MFEs)
print(results_df)

n1 = 7
n2 = 8

for seq_no in range (n1, n2):

    chain_file = seq_files[seq_no]
    chain_name_with_ext = os.path.basename(chain_file)        
    chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
    lp_file_name = chain_name_without_ext
    seq_data = parse_seq_file(chain_file)
    rna = seq_data['sequence']

    print(chain_name_without_ext)
    print(rna)
    print(len(rna))

    # subseq_len = round(len(this_RNA)*0.95)
    # s_start = len(this_RNA)-subseq_len+1

    # for i in range(s_start):
    #     gen_MFE_start, lp_name_start, opt_time_start = optimize_multi_start(seq_files, seq_no, lpstart_dir, solstart_dir, i, subseq_len)

    # gen_MFE_start, lp_name_start, opt_time_start = optimize_start(seq_files, seq_no, lpstart_dir, solstart_dir)

    # for i in range(1,s_start+1):
    #     f1_gen_start, fbeta_gen_start, MCC_gen_start, f1_rnastruct_start, fbeta_rnastruct_start, rna_len_start, MCC_rnastruct_start = sol_analyse(seq_files, seq_no, sol_dir, dot_bracket_start_dir, dot_bracket_archive_dir, dot_bracket_rnastructure_dir, 1, i)

    # for i in range(s_start):
    #     f1_gen, fbeta_gen, MCC_gen, f1_rnastruct, fbeta_rnastruct, MCC_rnastruct, f1_vienna, fbeta_vienna, MCC_vienna,rna_len = sol_analyse(seq_files, seq_no, sol_dir, dot_bracket_start_dir, dot_bracket_archive_dir, dot_bracket_rnastructure_dir, dot_bracket_viennaRNA_dir, 1, i)

    # gen_MFE, lp_name, opt_time = optimize_lilp(seq_files, seq_no, lp_dir, sol_dir, s_start, solstart_dir) 

    model_name = 'lilp'
    start_name = 'lilp-start'
    stem = True
    hairpin = True
    internal = True
    bulge = False
    multi = False
    start = True
    gen_MFE, lp_name, opt_time = optimize_lilp(rna, lp_file_name, model_name, stem, hairpin, internal, bulge, multi, lp_dir, incumbent_dir, sol_dir)
    # gen_MFE, lp_name, opt_time = optimize_lilp(rna, lp_file_name, model_name, stem, hairpin, internal, bulge, multi, lp_dir, incumbent_dir, sol_dir, start, start_name, solstart_dir)

    # f1_gen, fbeta_gen, MCC_gen, f1_rnastruct, fbeta_rnastruct, rna_len, MCC_rnastruct = sol_analyse(seq_files, seq_no, sol_dir,dot_bracket_dir, dot_bracket_archive_dir, dot_bracket_rnastructure_dir, 0)

    f1_gen, fbeta_gen, MCC_gen, f1_rnastruct, fbeta_rnastruct, MCC_rnastruct, f1_vienna, fbeta_vienna, MCC_vienna,rna_len = sol_analyse(seq_files, seq_no, sol_dir, model_name, dot_bracket_dir, dot_bracket_archive_dir, dot_bracket_rnastructure_dir, dot_bracket_viennaRNA_dir, 0)

    # write_results_to_file(lp_name, rna_len, opt_time, gen_MFE/100, ref_MFEs[seq_no], rna_MFEs[seq_no], round(f1_gen,2), round(f1_rnastruct,2), round(fbeta_gen,2), round(fbeta_rnastruct,2), round(MCC_gen,2), round(MCC_rnastruct,2), filename="ilp_LILP_0.2_results.txt")

    write_results_to_file(lp_name, rna_len, opt_time, gen_MFE/100, ref_MFEs[seq_no], rna_MFEs[seq_no], vienna_MFEs[seq_no], round(f1_gen,2), round(f1_rnastruct,2), round(f1_vienna,2), round(fbeta_gen,2), round(fbeta_rnastruct,2), round(fbeta_vienna,2), round(MCC_gen,2), round(MCC_rnastruct,2), round(MCC_vienna,2), filename="LILP_50_60_results.txt")
