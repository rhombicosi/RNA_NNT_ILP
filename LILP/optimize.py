import os
from pathlib import Path
import time
from utils.constants_paths import *
from utils.prepro_run import *
from lilp import *

def optimize(rna: str, model_name: str, stem: bool, hairpin: bool, internal: bool, bulge: bool, multi: bool, lp_dir: str, incumbent_dir: str, sol_dir: str, start = None, start_name = None, solstart_dir = None) -> None:
    
    rna_model = LILP(rna, model_name)
    rna_model.create_variables(stem, hairpin, internal, bulge, multi)
    rna_model.create_constraints(stem, hairpin, internal, bulge, multi)
    rna_model.create_objective(stem, hairpin, internal, bulge, multi)  
    
    rna_model.model.write(f'{lp_dir}/{lp_file_name}-{model_name}.lp')

    rna_model.model.setParam("LogFile", f'{grb_log_dir}/{lp_file_name}-log-{model_name}')
    rna_model.model.setParam(GRB.Param.SolFiles, f'{incumbent_dir}/{lp_file_name}-incumbent-{model_name}.lp')
    # rna_model.model.setParam("CliqueCuts", 1)
    # rna_model.model.setParam("RLTCuts", 1)
    # rna_model.model.setParam("ZeroHalfCuts", 1)
    # rna_model.model.setParam("RelaxLiftCuts", 2)
    # rna_model.model.setParam("MIPFocus", 2)
    rna_model.model.setParam("Heuristics", 0)

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

    if start:
        rna_model.model.NumStart = 1
        rna_model.model.update()
        solvars = read_sol(f'{solstart_dir}/{lp_file_name}-{start_name}.sol')

        # now set MIP start values using the Start attribute, e.g.:
        for v in rna_model.model.getVars(): 
            if v.VarName in solvars.keys():
                v.Start = round(int(solvars[v.VarName]), 1)
        rna_model.model.update()

    opt_start_time = time.time()
    rna_model.model.optimize()
    opt_time = time.time() - opt_start_time
    print(f'OPTIMIZATION TIME :: {opt_time}')

    rna_model.model.write(f'{sol_dir}/{lp_file_name}-{model_name}.sol')

    print(f'Obj: {rna_model.model.ObjVal:g}')

seq_number = 1

chain_file = seq_files[seq_number]
chain_name_with_ext = os.path.basename(chain_file)
chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
lp_file_name = chain_name_without_ext
seq_data = parse_seq_file(chain_file)

rna = seq_data['sequence']
print(chain_file)
print(rna)
print(len(rna))

start_name = 'lilp-start'
stem = True
hairpin = True
internal = False
bulge = False
multi = False
start = False
optimize(rna, start_name, stem, hairpin, internal, bulge, multi, lpstart_dir, incumbent_start_dir, solstart_dir)
model_name = 'lilp'
stem = True
hairpin = True
internal = True
bulge = True
multi = True
start = True
optimize(rna, model_name, stem, hairpin, internal, bulge, multi, lp_dir, incumbent_dir, sol_dir, start, start_name, solstart_dir)

# # process solution(s) to dot-bracket
# filepath = f'{sol_dir}/{lp_file_name}-lilp.sol'
# pairs2brackets(filepath, rna)

# calculate_sol_energy(filepath, rna)

# script_dir = os.path.dirname(os.path.abspath(__file__))
# parent_dir = os.path.dirname(script_dir)

# for f in range(12):
#     filepath = os.path.join(parent_dir, f'lilp_{seq_no}_incumbent__{f}.sol')
#     pairs2brackets(filepath, rna)

############################

# seq_number = 1

# chain_file = seq_files[seq_number]
# chain_name_with_ext = os.path.basename(chain_file)
# chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
# lp_file_name = chain_name_without_ext
# seq_data = parse_seq_file(chain_file)

# rna = seq_data['sequence']
# print(chain_file)
# print(rna)
# print(len(rna))

# start_model = LILP(rna, f'MIP-start-{seq_number}')
# start_model.create_variables(1, 1, 0, 0, 0)
# start_model.create_constraints(1, 1, 0, 0, 0)
# start_model.create_objective(1, 1, 0, 0, 0)

# start_model.model.write(f'{lpstart_dir}/{lp_file_name}-start.lp')

# start_model.model.setParam("LogFile", f'{grb_log_dir}/log-{lp_file_name}-start')
# # start_model.model.setParam("MIPFocus", 3)
# # start_model.model.setParam("Heuristics", 0)

# sorted_internals = sorted(start_model.internal_loops, key=lambda x: x.size)
# sorted_hairpins = sorted(start_model.hairpin_loops, key=lambda x: x.size)
# # sorted_bulges = sorted(start_model.bulge_loops, key=lambda x: x.size)
# # sorted_multis = sorted(start_model.multi_loops, key=lambda x: x.size)
# sorted_stems = sorted(start_model.stem_loops, key=lambda x: x.distance)

# for il in sorted_internals:
#     il.var.setAttr("BranchPriority", round(200/il.size))

# for hl in sorted_hairpins:
#     hl.var.setAttr("BranchPriority", round(200/hl.size))

# # for bl in sorted_bulges:
# #     bl.var.setAttr("BranchPriority", round(100/hl.size))

# # for ml in sorted_multis:
# #     ml.var.setAttr("BranchPriority", round(100/hl.size))

# for sl in sorted_stems:
#     sl.var.setAttr("BranchPriority", 10 * sl.distance)

# start_model.model.setParam(GRB.Param.SolFiles, f'{incumbent_start_dir}/{lp_file_name}-lilp-incumbent-start.lp')

# opt_start_time = time.time()
# start_model.model.optimize()
# opt_time = time.time() - opt_start_time

# start_model.model.write(f'{solstart_dir}/{lp_file_name}-lilp-start.sol')

# print(f'Obj: {start_model.model.ObjVal:g}')

# ##################


# rna_model = LILP(rna, f'MIP-{seq_number}')
# rna_model.create_variables(1, 1, 1, 1, 1)
# rna_model.create_constraints(1, 1, 1, 1, 1)
# rna_model.create_objective(1, 1, 1, 1, 1)

# rna_model.model.write(f'{lp_dir}/{lp_file_name}-lilp.lp')

# rna_model.model.setParam("LogFile", f'{grb_log_dir}/{lp_file_name}-log-lilp')
# # rna_model.model.setParam("CliqueCuts", 1)
# # rna_model.model.setParam("RLTCuts", 1)
# # rna_model.model.setParam("ZeroHalfCuts", 1)
# # rna_model.model.setParam("RelaxLiftCuts", 2)
# # rna_model.model.setParam("MIPFocus", 2)
# rna_model.model.setParam("Heuristics", 0)

# sorted_internals = sorted(rna_model.internal_loops, key=lambda x: x.size)
# sorted_hairpins = sorted(rna_model.hairpin_loops, key=lambda x: x.size)
# sorted_bulges = sorted(rna_model.bulge_loops, key=lambda x: x.size)
# sorted_multis = sorted(rna_model.multi_loops, key=lambda x: x.size)
# sorted_stems = sorted(rna_model.stem_loops, key=lambda x: x.distance)

# for il in sorted_internals:
#     il.var.setAttr("BranchPriority", round(200/il.size))

# for hl in sorted_hairpins:
#     hl.var.setAttr("BranchPriority", round(200/hl.size))

# for bl in sorted_bulges:
#     bl.var.setAttr("BranchPriority", round(100/hl.size))

# for ml in sorted_multis:
#     ml.var.setAttr("BranchPriority", round(100/hl.size))

# for sl in sorted_stems:
#     sl.var.setAttr("BranchPriority", 10 * sl.distance)

# rna_model.model.setParam(GRB.Param.SolFiles, f'{incumbent_dir}/{lp_file_name}-lilp-incumbent.lp')

# rna_model.model.NumStart = 1
# rna_model.model.update()
# solvars = read_sol(f'{solstart_dir}/{lp_file_name}-lilp-start.sol')

# # now set MIP start values using the Start attribute, e.g.:
# for v in rna_model.model.getVars(): 
#     if v.VarName in solvars.keys():
#         v.Start = round(int(solvars[v.VarName]), 1)
# rna_model.model.update()

# opt_start_time = time.time()
# rna_model.model.optimize()
# opt_time = time.time() - opt_start_time

# # rna_model.model.write(f'lilp-{seq_number}.sol')
# rna_model.model.write(f'{sol_dir}/{lp_file_name}-lilp.sol')

# print(f'Obj: {rna_model.model.ObjVal:g}')






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

# counter = 0

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

