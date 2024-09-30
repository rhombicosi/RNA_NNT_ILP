import pandas as pd
from sol_converter import *
from prepro_run import *
from constants_paths import *
from ilp_grb import *

def add_column(df, column_name, values):
    df[column_name] = values
    return df

ref_mfe_files = get_filenames(efn2_archive_dir, '.txt')
ref_mfe_names = [os.path.splitext(os.path.basename(f))[0] for f in ref_mfe_files]
results_df = pd.DataFrame(index=ref_mfe_names)

ref_MFEs = [get_energy_from_ct_file(f) for f in ref_mfe_files]

rna_mfe_names = get_filenames(ct_rnastruct_dir, '.ct')
rna_MFEs = [get_energy_from_ct_file(f) for f in rna_mfe_names]

add_column(results_df, 'MFE_ref', ref_MFEs)
add_column(results_df, 'MFE_rna', rna_MFEs)
print(results_df)

n1 = 47
n2 = 48

for seq_no in range (n1, n2):

    gen_MFE, lp_name = optimize(seq_files, seq_no, lp_dir, sol_dir)

    f1_ref, fbeta_ref, f1_rnastruct, fbeta_rnastruct, rna_len = sol_analyse(seq_files, seq_no, sol_dir, dot_bracket_archive_dir, dot_bracket_rnastructure_dir)

    write_results_to_file(lp_name, rna_len, gen_MFE/100, ref_MFEs[seq_no], rna_MFEs[seq_no], round(f1_ref,2), round(f1_rnastruct,2), round(fbeta_ref,2), round(fbeta_rnastruct,2), filename="ilp_multi_results.txt")
