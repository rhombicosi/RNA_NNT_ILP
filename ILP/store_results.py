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

rna_mfe_files = get_filenames(ct_rnastruct_dir, '.ct')
rna_MFEs = [get_energy_from_ct_file(f) for f in rna_mfe_files]

add_column(results_df, 'MFE_ref', ref_MFEs)
add_column(results_df, 'MFE_rna', rna_MFEs)
print(results_df)

for seq in range (9, 10):
    obj, lp_name = optimize(seq_files, seq, lp_dir, sol_dir)

    lp_files.append(lp_name)
    MFE_gen.append(obj)

print(lp_files)
print(MFE_gen)

for seq in range (9, 10):
    f1_ref_val, f1beta_ref_val, f1_rnastruct_val, f1beta_rnastruct_val = sol_analyse(seq_files, seq, sol_dir, dot_bracket_archive_dir, dot_bracket_rnastructure_dir)

    f1_ref.append(f1_ref_val)
    f1beta_ref.append(f1beta_ref_val)
    f1_rnastruct.append(f1_rnastruct_val)
    f1beta_rnastruct.append(f1beta_rnastruct_val)

print(f1_ref)
print(f1beta_ref)
print(f1_rnastruct)
print(f1beta_rnastruct)

print(results_df.iloc[9])
