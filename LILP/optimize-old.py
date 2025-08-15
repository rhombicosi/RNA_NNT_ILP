import os
from pathlib import Path
import time
from utils.constants_paths import *
from utils.prepro_run import *
from lilp import *


seq_len = 60
seq_no = 0

chain_file = seq_files[seq_no]
chain_name_with_ext = os.path.basename(chain_file)        
chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
lp_file_name = chain_name_without_ext
seq_data = parse_seq_file(chain_file)

rna = seq_data['sequence']
print(rna)

rna_start_model = LILP(rna, 'lilp-test')
rna_start_model.create_base_pairs()
rna_start_model.create_stem_loops()
rna_start_model.create_hairpin_loops()
rna_start_model.create_internal_loops()
# rna_start_model.create_bulge_loops()
# rna_start_model.create_multi_loops()

rna_start_model.add_single_pair_constraints()
rna_start_model.add_no_crossing_constraints()

rna_start_model.add_stem_constraints()

rna_start_model.create_nucleotides()
rna_start_model.add_unpaired_nucleotides_constraints()

rna_start_model.add_hairpin_size_constraints()
rna_start_model.add_hairpin_ifthen_constraints()
# rna_model.add_hairpin_onlyif_constraints()
# rna_model.add_hairpin_max_number_constraint()

rna_start_model.add_internal_size_constraints()
rna_start_model.add_internal_ifthen_constraints()
rna_start_model.add_internal_onlyif_constraints()
# rna_model.add_internal_max_number_constraint()
# rna_model.model.addConstr(rna_model.model.getVarByName(f'INTERNAL_5_39_11_36') == 1)

# rna_start_model.add_bulge_size_constraints()
# rna_start_model.add_bulge_ifthen_constraints()
# # rna_model.add_bulge_onlyif_constraints()
# # rna_model.add_bulge_max_number_constraint()

# rna_start_model.add_multi_size_constraints()
# rna_start_model.add_multi_ifthen_constraints()
# # rna_model.add_multi_onlyif_constraints()
# # rna_model.add_multi_max_number_constraint()

rna_start_model.create_objective(1, 1, 1, 0, 0)

rna_start_model.model.write(f'lilp-start-{seq_no}.lp')

for il in rna_start_model.internal_loops:
    il.var.setAttr("BranchPriority", round(200/il.size))

for hl in rna_start_model.hairpin_loops:
    hl.var.setAttr("BranchPriority", round(200/hl.size))

# for bl in rna_start_model.bulge_loops:
#     bl.var.setAttr("BranchPriority", round(100/hl.size))

# for ml in rna_start_model.multi_loops:
#     ml.var.setAttr("BranchPriority", round(100/hl.size))

# for sl in rna_start_model.stem_loops:
#     sl.var.setAttr("BranchPriority", 10 * sl.distance)

opt_start_time = time.time()
rna_start_model.model.optimize()
opt_time = time.time() - opt_start_time

rna_start_model.model.write(f'lilp-start-{seq_no}.sol')

print(f'Obj: {rna_start_model.model.ObjVal:g}')

script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
filepath = os.path.join(parent_dir, f'lilp-start-{seq_no}.sol')

pairs2brackets(filepath, rna)
# calculate_sol_energy(filepath, rna)

######################

rna_model = LILP(rna, 'lilp-test')
rna_model.create_base_pairs()
rna_model.add_single_pair_constraints()
rna_model.add_no_crossing_constraints()

rna_model.create_stem_loops()
rna_model.add_stem_constraints()

rna_model.create_hairpin_loops()
rna_model.create_internal_loops()
rna_model.create_bulge_loops()
rna_model.create_multi_loops()

rna_model.create_nucleotides()
rna_model.add_unpaired_nucleotides_constraints()

rna_model.add_hairpin_size_constraints()
rna_model.add_hairpin_ifthen_constraints()
# rna_model.add_hairpin_onlyif_constraints()
# rna_model.add_hairpin_max_number_constraint()

rna_model.add_internal_size_constraints()
rna_model.add_internal_ifthen_constraints()
rna_model.add_internal_onlyif_constraints()
# rna_model.add_internal_max_number_constraint()
# rna_model.model.addConstr(rna_model.model.getVarByName(f'INTERNAL_5_39_11_36') == 1)

rna_model.add_bulge_size_constraints()
rna_model.add_bulge_ifthen_constraints()
# rna_model.add_bulge_onlyif_constraints()
# rna_model.add_bulge_max_number_constraint()

rna_model.add_multi_size_constraints()
rna_model.add_multi_ifthen_constraints()
# rna_model.add_multi_onlyif_constraints()
# rna_model.add_multi_max_number_constraint()

rna_model.create_objective(1, 1, 1, 1, 1)

rna_model.model.write(f'lilp-{seq_no}.lp')

for il in rna_model.internal_loops:
    il.var.setAttr("BranchPriority", round(200/il.size))

for hl in rna_model.hairpin_loops:
    hl.var.setAttr("BranchPriority", round(200/hl.size))

# for bl in rna_model.bulge_loops:
#     bl.var.setAttr("BranchPriority", round(100/hl.size))

# for ml in rna_model.multi_loops:
#     ml.var.setAttr("BranchPriority", round(100/hl.size))

# for sl in rna_model.stem_loops:
#     sl.var.setAttr("BranchPriority", 10 * sl.distance)


rna_model.model.NumStart = 1
rna_model.model.update()
solvars = read_sol(f'lilp-start-{seq_no}.sol')

# start values
for v in rna_model.model.getVars(): 
    if v.VarName in solvars.keys():
        v.Start = round(int(solvars[v.VarName]), 1)
rna_model.model.update()

opt_start_time = time.time()
rna_model.model.optimize()
opt_time = time.time() - opt_start_time

rna_model.model.write(f'lilp-{seq_no}.sol')

print(f'Obj: {rna_model.model.ObjVal:g}')


# script_dir = os.path.dirname(os.path.abspath(__file__))
# parent_dir = os.path.dirname(script_dir)
# filepath = os.path.join(parent_dir, f'lilp-{seq_no}.sol')


# pairs2brackets(filepath, rna)

# calculate_sol_energy(filepath, rna)