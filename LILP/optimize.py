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

        # start values
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

# start_name = 'lilp-start'
# stem = True
# hairpin = True
# internal = False
# bulge = False
# multi = False
# start = False
# optimize(rna, start_name, stem, hairpin, internal, bulge, multi, lpstart_dir, incumbent_start_dir, solstart_dir)
model_name = 'lilp'
stem = True
hairpin = True
internal = True
bulge = True
multi = True
# start = True
optimize(rna, model_name, stem, hairpin, internal, bulge, multi, lp_dir, incumbent_dir, sol_dir)
# optimize(rna, model_name, stem, hairpin, internal, bulge, multi, lp_dir, incumbent_dir, sol_dir, start, start_name, solstart_dir)

# process solution(s) to dot-bracket
filepath = f'{sol_dir}/{lp_file_name}-{model_name}.sol'
pairs2brackets(filepath, rna)

# calculate_sol_energy(filepath, rna)

# script_dir = os.path.dirname(os.path.abspath(__file__))
# parent_dir = os.path.dirname(script_dir)

# for f in range(12):
#     filepath = os.path.join(parent_dir, f'lilp_{seq_no}_incumbent__{f}.sol')
#     pairs2brackets(filepath, rna)
