from lilp import *


seq_len = 60
seq_no = 1

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

rna_model = LILP(rna)
rna_model.create_variables(1, 1, 1, 1, 1)
rna_model.create_constraints(1, 1, 1, 1, 1)
rna_model.create_objective(1, 1, 1, 1, 1)

rna_model.model.write(f'lilp-{seq_no}.lp')

rna_model.model.setParam("LogFile", f'log-{seq_no}')
# rna_model.model.setParam("MIPFocus", 2)
# rna_model.model.setParam("Heuristics", 0)

sorted_internals = sorted(rna_model.internal_loops, key=lambda x: x.size)
sorted_hairpins = sorted(rna_model.hairpin_loops, key=lambda x: x.size)
sorted_bulges = sorted(rna_model.bulge_loops, key=lambda x: x.size)
sorted_multis = sorted(rna_model.multi_loops, key=lambda x: x.size)
sorted_stems = sorted(rna_model.stem_loops, key=lambda x: x.distance)

for il in sorted_internals:
    il.var.setAttr("BranchPriority", round(200/il.size))

for hl in sorted_hairpins:
    hl.var.setAttr("BranchPriority", round(200/hl.size))

for bl in sorted_bulges:
    bl.var.setAttr("BranchPriority", round(100/hl.size))

for ml in sorted_multis:
    ml.var.setAttr("BranchPriority", round(100/hl.size))

for sl in sorted_stems:
    sl.var.setAttr("BranchPriority", 10 * sl.distance)

rna_model.model.setParam(GRB.Param.SolFiles, f'lilp_{seq_no}_incumbent')

opt_start_time = time.time()
rna_model.model.optimize()
opt_time = time.time() - opt_start_time

rna_model.model.write(f'lilp-{seq_no}.sol')

print(f'Obj: {rna_model.model.ObjVal:g}')


script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
filepath = os.path.join(parent_dir, f'lilp-{seq_no}.sol')

pairs2brackets(filepath, rna)

# calculate_sol_energy(filepath, rna)

# script_dir = os.path.dirname(os.path.abspath(__file__))
# parent_dir = os.path.dirname(script_dir)

# for f in range(12):
#     filepath = os.path.join(parent_dir, f'lilp_{seq_no}_incumbent__{f}.sol')
#     pairs2brackets(filepath, rna)





# rna = "GCCGCGAACCCCGCCAGGCCCGGAAGGGAGCAACGGUAGUGGUGGAU"
# bp1 = BasePair(1,43,rna)
# bp2 = BasePair(2,42,rna)
# bp3 = BasePair(3,41,rna)
# bp4 = BasePair(4,40,rna)
# bp5 = BasePair(5,39,rna)
# bp6 = BasePair(11,36,rna)
# bp7 = BasePair(12,35,rna)
# bp8 = BasePair(13,34,rna)
# bp9 = BasePair(19,28,rna)
# bp10 = BasePair(20,27,rna)
# bp11 = BasePair(21,26,rna)
# loop1 = StemLoop((bp1,bp2),rna)
# loop2 = StemLoop((bp2,bp3),rna)
# loop3 = StemLoop((bp3,bp4),rna)
# loop4 = StemLoop((bp4,bp5),rna)
# loop5 = InternalLoop((bp5, bp6), rna)
# loop6 = StemLoop((bp6,bp7),rna)
# loop7 = StemLoop((bp7,bp8),rna)
# loop8 = InternalLoop((bp8, bp9), rna)
# loop9 = StemLoop((bp9,bp10),rna)
# loop10 = StemLoop((bp10,bp11),rna)
# loop11 = HairpinLoop([bp11],rna)

# print(f'{loop1.energy}')
# print(f'{loop2.energy}')
# print(f'{loop3.energy}')
# print(f'{loop4.energy}')
# print(f'{loop5.energy}')
# print(f'{loop6.energy}')
# print(f'{loop7.energy}')
# print(f'{loop8.energy}')
# print(f'{loop9.energy}')
# print(f'{loop10.energy}')
# print(f'{loop11.energy}')

counter = 0

# print(type(rna_model.base_pairs))

# for bp in rna_model.base_pairs:    
#     print(f"[{counter}]::Base pair: ({bp.i}, {bp.j}), First pair {bp.nt1 + bp.nt2} energy: {bp.pair_penalty_energy}, Variable name: {bp.var.VarName}, Gurobi var: {bp.var}")
#     counter+=1

# for bp in rna_model.first_pairs:    
#     print(f"[{counter}]::Base pair: ({bp.i}, {bp.j}), Variable name: {bp.var.VarName}, Gurobi var: {bp.var}")
#     counter+=1

# for h in rna_model.hairpin_loops:    
#     print(f"[{counter}]:: Variable name: {h.var.VarName}, Energy: {h.energy}, Gurobi var: {h.var}")
#     counter+=1

# for s in rna_model.stem_loops:    
#     print(f"[{counter}]:: Variable name: {s.var.VarName}, Energy of {rna[s.first_pair.i - 1] + rna[s.first_pair.j - 1]} followed by {rna[s.last_pair.i - 1] + rna[s.last_pair.j - 1]}: {s.energy}, Gurobi var: {s.var}")
#     counter+=1

# for i in rna_model.internal_loops:
#     if i.is_valid_size():  
#         print(f"[{counter}]:: Variable name: {i.var.VarName}, Size: {i.size}, Energy {i.subtype}: {i.energy}, Gurobi var: {i.var}")
#     counter+=1

# loop = rna_model.internal_loops[93]
# print(loop.var.VarName)

# for b in rna_model.bulge_loops:    
#     print(f"[{counter}]:: Variable name: {b.var.VarName}, Energy: {b.energy}, Gurobi var: {b.var}")
#     counter+=1

# for m in rna_model.multi_loops:    
#     print(f"[{counter}]:: Variable name: {m.var.VarName}, Size: {m.size}, Energy: {m.energy}, Gurobi var: {m.var}")
#     counter+=1

# dloop = Loop([rna_model.base_pairs[67], rna_model.base_pairs[90]],rna)

