from lilp import *


seq_len = 60
seq_no = 0

cwd = Path.cwd()
code_path = Path(__file__).parent.parent
arch_rel_path = '../../ARCHIVE II/'
archive_path = (code_path/arch_rel_path).resolve()
seq_len_dir = f'RNA_seq_{seq_len}'

chain_dir = os.path.join(archive_path, seq_len_dir)
seq_files = get_filenames(chain_dir, '.seq')

chain_file = seq_files[seq_no]
chain_name_with_ext = os.path.basename(chain_file)        
chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
lp_file_name = chain_name_without_ext
seq_data = parse_seq_file(chain_file)

rna = seq_data['sequence']
print(rna)

rna_model = LILPModel(rna)
rna_model.create_base_pairs()
rna_model.create_stem_loops()
rna_model.create_hairpin_loops()
rna_model.create_internal_loops()
rna_model.create_bulge_loops()
rna_model.create_multi_loops()

rna_model.add_single_pair_constraints()
rna_model.add_no_crossing_constraints()

rna_model.add_stem_constraints()

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


var_STEM = [var for var in rna_model.model.getVars() if 'STEM' in var.VarName]
for i in range(len(var_STEM)):
    var_STEM[i].setAttr("BranchPriority",100)

opt_start_time = time.time()
rna_model.model.optimize()
opt_time = time.time() - opt_start_time

rna_model.model.write(f'lilp-{seq_no}.sol')

print(f'Obj: {rna_model.model.ObjVal:g}')


script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
filepath = os.path.join(parent_dir, f'lilp-{seq_no}.sol')

pairs2brackets(filepath, rna)

calculate_sol_energy(filepath, rna)