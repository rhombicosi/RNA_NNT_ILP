##### needed only to get sequence for test #####
import os
from pathlib import Path
import pandas as pd
from prepro_utils import *
import gurobipy as gp
from gurobipy import GRB
from lilp_config import *
from basepair import *
from dloop import *
from stemloop import *

class LILPModel:
    def __init__(self, rna_seq: str):
        self.rna_seq = rna_seq
        self.model = gp.Model(f'MIP')
        self.base_pairs : List[BasePair] = []
        self.nucleotides = []
        self.first_pairs = []
        self.last_pairs = []
        self.hairpin_loops = []
        self.stem_loops = []
        self.internal_loops = []
        self.bulge_loops = []
        self.multi_loops = []    

    def create_base_pairs(self) -> None:
        n = len(self.rna_seq)

        for i in range(1, n): 
            for j in range(i + MIN_D + 1, n + 1):
                bp = BasePair(i, j, self.rna_seq)
                var = bp.add_variable(self.model, 'P')
                if var:
                    self.base_pairs.append(bp)
        self.model.update()

    def create_first_pairs(self) -> None:

        for bp in self.base_pairs:
            fp = BasePair(bp.i, bp.j, self.rna_seq)
            self.first_pairs.append(fp)
            fp.add_variable(self.model, 'F') 
        self.model.update()

    def create_last_pairs(self) -> None:

        for bp in self.base_pairs:
            lp = BasePair(bp.i, bp.j, self.rna_seq)
            self.last_pairs.append(lp)
            lp.add_variable(self.model,'L')
        self.model.update()

    def create_nucleotides(self) -> None:
        n = len(self.rna_seq)

        for i in range(1, n + 1):
            var = self.model.addVar(vtype=GRB.BINARY, name=f'X_{i}')   
            self.nucleotides.append(var)     
        self.model.update()

    def create_hairpin_loops(self) -> None:

        for bp in self.base_pairs:
            hairpin = Loop([bp], self.rna_seq)
            hairpin.add_variable(self.model)
            self.hairpin_loops.append(hairpin)
        self.model.update() 

    def create_stem_loops(self) -> None:

        for bp1 in self.base_pairs:
            bp2 = BasePair._find_base_pairs_matches(self.base_pairs, bp1.i + 1, bp1.j - 1)
            if bp2:
                stem = StemLoop([bp1, bp2], self.rna_seq)
                stem.add_variable(self.model)
                self.stem_loops.append(stem)
        self.model.update()

    def create_internal_loops(self) -> None:

        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                if bp2.i > bp1.i + 1 and bp2.j < bp1.j - 1:
                    internal = Loop([bp1, bp2], self.rna_seq)
                    internal.add_variable(self.model)
                    self.internal_loops.append(internal)
        self.model.update()

    def create_bulge_loops(self) -> None:

        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                if (bp2.i == bp1.i + 1 and bp2.j < bp1.j - 1) or (bp2.i > bp1.i + 1 and bp2.j == bp1.j - 1):
                    bulge = Loop([bp1, bp2], self.rna_seq)
                    bulge.add_variable(self.model)
                    self.bulge_loops.append(bulge)
        self.model.update()

    def create_multi_loops(self) -> None:

        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                for bp3 in self.base_pairs:
                    if bp2.i > bp1.i and bp3.i > bp2.j and bp1.j > bp3.j:
                        multi = Loop([bp1, bp2, bp3], self.rna_seq)
                        multi.add_variable(self.model)
                        self.multi_loops.append(multi)
        self.model.update()

    def add_single_pair_constraints(self) -> None:
        n = len(self.rna_seq)

        for i in range(1, n + 1):
            BasePair.create_single_pair_constraint(self.model, self.base_pairs, i)
        self.model.update()

    def add_no_crossing_constraints(self) -> None: 
        for bp1 in self.base_pairs: 
            for bp2 in self.base_pairs:
                BasePair.create_no_crossing_constraint(self.model, bp1, bp2)
        self.model.update()

    # def create_no_crossing_constraints(self) -> None: 

    #     for bp1 in self.base_pairs: 
    #         for bp2 in self.base_pairs:
    #             if bp2.i > bp1.i and bp2.i < bp1.j and bp2.j > bp1.j:
    #                 inequality = gp.LinExpr(0)
    #                 inequality.add(gp.LinExpr([1.0,1.0],[bp1.var,bp2. var]))
    #                 self.model.addConstr(inequality <= 1, f'NC-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}')
    #     self.model.update()

    def add_stem_constraints(self) -> None:
        for sl in self.stem_loops:
            sl.create_stem_constraints(self.model)
        self.model.update()
        
    def add_first_pair_constraints(self) -> None:
        for sl in self.stem_loops:
            sl.create_first_pair_constraints(self.model, self.stem_loops, self.first_pairs)
        self.model.update()

    def add_last_pair_constraints(self) -> None:
        for sl in self.stem_loops:    
            sl.create_last_pair_constraints(self.model, self.stem_loops, self.last_pairs)
        self.model.update()

    def create_unpaired_nucleotides_constraints(self) -> None:
        n = len(self.rna_seq)

        for i in range(1, n + 1):
            inequality = gp.LinExpr(0)
            matches = BasePair._find_base_pairs_with_index(self.base_pairs, i)

            if matches:
                for bp in matches:
                    inequality.add(gp.LinExpr([1.0],[bp.var]))
                
                inequality.add(gp.LinExpr([1.0],[self.nucleotides[i-1]])) 
                self.model.addConstr(inequality == 1, f'UN-{i}')
        self.model.update()

    def create_hairpin_size_constraints(self) -> None:

        for hl in self.hairpin_loops:
            if not hl.is_valid_size():
                inequality = gp.LinExpr([1],[hl.var])
                self.model.addConstr(inequality == 0, f'HS-{hl.base_pairs[0].i}-{hl.base_pairs[0].j}')
        self.model.update()

    def create_hairpin_ifthen_constraints(self) -> None:

        for hl in self.hairpin_loops:
            inequality = gp.LinExpr(0)
            for u in range(hl.base_pairs[0].i + 1, hl.base_pairs[0].j):
                inequality.add(gp.LinExpr([1],[self.nucleotides[u-1]]))
            
            inequality.add(gp.LinExpr([1, -1],[hl.base_pairs[0].var, hl.var]))            
            self.model.addConstr(inequality <= hl.size, f'HIT-{hl.base_pairs[0].i}-{hl.base_pairs[0].j}')
        self.model.update()

    def create_hairpin_onlyif_constraints(self) -> None:

        for hl in self.hairpin_loops:
            for u in range(hl.base_pairs[0].i + 1, hl.base_pairs[0].j):
                inequality = gp.LinExpr([2],[hl.var])
                matches = self._find_base_pairs_with_index(self.base_pairs, u)

                for bp in matches:
                    inequality.add(gp.LinExpr([1],[bp.var]))
                
                inequality.add(gp.LinExpr([-1], [hl.base_pairs[0].var]))
                self.model.addConstr(inequality <= 1, f'HOI-{hl.base_pairs[0].i}-{hl.base_pairs[0].j}-{u}')
        self.model.update()

    def create_hairpin_max_number_constraints(self) -> None:
        inequality = gp.LinExpr(0)

        for hl in self.hairpin_loops:
            inequality.add(gp.LinExpr([1],[hl.var]))
        self.model.addConstr(inequality <= MAX_NUM_OF_LOOPS[hl.type], f'HMN')
        self.model.update()

    def create_internal_size_constraints(self) -> None:

        for il in self.internal_loops:
            if not il.is_valid_size():
                inequality = gp.LinExpr([1],[il.var])
                self.model.addConstr(inequality == 0, f'IS-{il.base_pairs[0].i}-{il.base_pairs[0].j}-{il.base_pairs[1].i}-{il.base_pairs[1].j}')
        self.model.update()

#     def create_internal_ifthen_constaints(self) -> None:



# def internalIfThenConstraints(RNA, mip):
#     n = len(RNA)
#     for i in range(1, n - minI - 1 - minD - 1 - minI - 1):
#         for k in range(i + minI + 1, n - minI - 1 - minD - 1):
#             for l in range(k + minD + 1, n - minI  - 1):
#                 for j in range(l + minI + 1, n + 1):
#                     if RNA[i-1] + RNA[j-1] in cbp_list and RNA[k-1] + RNA[l-1] in cbp_list:
#                         inequality = gp.LinExpr(0)                        

#                         for u in range(i+1,k):
#                             inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))

#                         for u in range(l+1,j):
#                             inequality.add(gp.LinExpr([1],[mip.getVarByName(f'X({u})')]))
                            
#                         inequality.add(gp.LinExpr([1,1,-1],[mip.getVarByName(f'P({k},{l})'),mip.getVarByName(f'P({i},{j})'),mip.getVarByName(f'I({i},{k},{l},{j})')]))

#                         mip.addConstr(inequality <= k-i+j-l-1, f'CIIFT{i}-{k}-{l}-{j}')

#     return inequality


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
rna_model = LILPModel(rna)
rna_model.create_base_pairs()
rna_model.create_first_pairs()
rna_model.create_last_pairs()
rna_model.create_hairpin_loops()
rna_model.create_stem_loops()
# rna_model.create_internal_loops()
# rna_model.create_bulge_loops()
# rna_model.create_multi_loops()
rna_model.add_single_pair_constraints()
rna_model.add_no_crossing_constraints()
rna_model.add_stem_constraints()
rna_model.add_first_pair_constraints()
rna_model.add_last_pair_constraints()
rna_model.create_nucleotides()
rna_model.create_unpaired_nucleotides_constraints()
# rna_model.create_hairpin_size_constraints()
# rna_model.create_hairpin_ifthen_constraints()
# rna_model.create_hairpin_onlyif_constraints()
# rna_model.create_hairpin_max_number_constraints()
# rna_model.create_internal_size_constraints()

rna_model.model.write(f'lilp-{seq_no}.lp')

# counter = 0

# print(type(rna_model.base_pairs))

# for bp in rna_model.base_pairs:    
#     print(f"[{counter}]::Base pair: ({bp.i}, {bp.j}), Variable name: {bp.var.VarName}, Gurobi var: {bp.var}")
#     counter+=1

# for bp in rna_model.first_pairs:    
#     print(f"[{counter}]::Base pair: ({bp.i}, {bp.j}), Variable name: {bp.var.VarName}, Gurobi var: {bp.var}")
#     counter+=1

# for h in rna_model.hairpin_loops:    
#     print(f"[{counter}]:: Variable name: {h.var.VarName}, Gurobi var: {h.var}")
#     counter+=1

# for s in rna_model.stem_loops:    
#     print(f"[{counter}]:: Variable name: {s.var.VarName}, Gurobi var: {s.var}")
#     counter+=1

# for i in rna_model.internal_loops:    
#     print(f"[{counter}]:: Variable name: {i.var.VarName}, Gurobi var: {i.var}")
#     counter+=1

# for b in rna_model.bulge_loops:    
#     print(f"[{counter}]:: Variable name: {b.var.VarName}, Gurobi var: {b.var}")
#     counter+=1

# for m in rna_model.multi_loops:    
#     print(f"[{counter}]:: Variable name: {m.var.VarName}, Gurobi var: {m.var}")
#     counter+=1

dloop = Loop([rna_model.base_pairs[67], rna_model.base_pairs[90]],rna)

# P :: 113 :::: Q :: 24 :::: H :: 113 :::: I :: 1053 :::: B :: 387 :::: M :: 2694

