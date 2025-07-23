##### needed only to get sequence for test #####
import os
from pathlib import Path
import pandas as pd
from prepro_utils import *
import gurobipy as gp
from gurobipy import GRB
from lilp_config import *
from basepair import *
from dloops import *


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
        
    def _find_base_pairs_with_index(self, base_pairs: List[BasePair], index: int) -> List[BasePair]:
        return [bp for bp in base_pairs if bp.i == index or bp.j == index]
    
    def _find_base_pairs_matches(self, base_pair: BasePair) -> BasePair:
        return next((bp for bp in self.base_pairs if bp.i == base_pair.i + 1 and bp.j == base_pair.j - 1), None)
    
    def _find_stem_matches(self, stem: RNALoop) -> RNALoop:
        bp = stem.base_pairs[0]
        return next((sl for sl in self.stem_loops if sl.base_pairs[1].i == bp.i and sl.base_pairs[1].j == bp.j), None)

    def create_base_pairs(self) -> None:
        n = len(self.rna_seq)

        for i in range(1, n): 
            for j in range(i + MIN_D + 1, n + 1):
                bp = BasePair(i, j, self.rna_seq)
                var = bp.add_variable(self.model, 'P')
                if var:
                    self.base_pairs.append(bp)
        self.model.update()

    def create_fp_pairs(self) -> None:

        for bp in self.base_pairs:
            fp = BasePair(bp.i, bp.j, self.rna_seq)
            self.last_pairs.append(bp)
            fp.add_variable(self.model) 
        self.model.update()

    def create_last_pairs(self) -> None:

        for bp in self.base_pairs:
            lp = BasePair(bp.i, bp.j, self.rna_seq)
            self.last_pairs.append(bp)
            lp.add_variable(self.model)
        self.model.update()

    def create_nucleotides(self) -> None:
        n = len(self.rna_seq)

        for i in range(1, n + 1):
            var = self.model.addVar(vtype=GRB.BINARY, name=f'X_{i}')   
            self.nucleotides.append(var)     
        self.model.update()

    def create_hairpin_loops(self) -> None:

        for bp in self.base_pairs:
            hairpin = RNALoop([bp], self.rna_seq)
            hairpin.add_variable(self.model)
            self.hairpin_loops.append(hairpin)
        self.model.update() 

    def create_stem_loops(self) -> None:

        for bp1 in self.base_pairs:
            bp2 = self._find_base_pairs_matches(bp1)
            if bp2:
                stem = RNALoop([bp1, bp2], self.rna_seq)
                stem.add_variable(self.model)
                self.stem_loops.append(stem)
        self.model.update()

    def create_internal_loops(self) -> None:

        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                if bp2.i > bp1.i + 1 and bp2.j < bp1.j - 1:
                    internal = RNALoop([bp1, bp2], self.rna_seq)
                    internal.add_variable(self.model)
                    self.internal_loops.append(internal)
        self.model.update()

    def create_bulge_loops(self) -> None:

        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                if (bp2.i == bp1.i + 1 and bp2.j < bp1.j - 1) or (bp2.i > bp1.i + 1 and bp2.j == bp1.j - 1):
                    bulge = RNALoop([bp1, bp2], self.rna_seq)
                    bulge.add_variable(self.model)
                    self.bulge_loops.append(bulge)
        self.model.update()

    def create_multi_loops(self) -> None:

        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                for bp3 in self.base_pairs:
                    if bp2.i > bp1.i and bp3.i > bp2.j and bp1.j > bp3.j:
                        multi = RNALoop([bp1, bp2, bp3], self.rna_seq)
                        multi.add_variable(self.model)
                        self.multi_loops.append(multi)
        self.model.update()

    def create_single_pair_constraints(self) -> None:
        n = len(self.rna_seq)

        for i in range(1, n + 1):
            inequality = gp.LinExpr(0)
            matches = self._find_base_pairs_with_index(self.base_pairs, i)

            if matches:
                for bp in matches:
                    inequality.add(gp.LinExpr([1.0], [bp.var]))
                self.model.addConstr(inequality <= 1, f'SP-{i}')
        self.model.update()

    def create_no_crossing_constraints(self) -> None: 

        for bp1 in self.base_pairs: 
            for bp2 in self.base_pairs:
                if bp2.i > bp1.i and bp2.i < bp1.j and bp2.j > bp1.j:
                    inequality = gp.LinExpr(0)
                    inequality.add(gp.LinExpr([1.0,1.0],[bp1.var,bp2. var]))
                    self.model.addConstr(inequality <= 1, f'NC-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}')
        self.model.update()

    def create_stem_constraints(self) -> None:

        for sl in self.stem_loops:
            ifthen_inequality = gp.LinExpr([2,-1,-1],[sl.var, sl.base_pairs[0].var, sl.base_pairs[1].var])
            self.model.addConstr(ifthen_inequality <= 0, f'SLIT-{sl.base_pairs[0].i}-{sl.base_pairs[0].j}')

            onlyif_inequality = gp.LinExpr([1,1,-1],[sl.base_pairs[0].var, sl.base_pairs[1].var, sl.var])
            self.model.addConstr(onlyif_inequality <= 1, f'SLOI-{sl.base_pairs[0].i}-{sl.base_pairs[0].j}')
        self.model.update()
        
    def create_first_pair_constraints(self) -> None:
        n = len(self.rna_seq)
        for sl1 in self.stem_loops:
            if sl1.base_pairs[0].i > 1 and sl1.base_pairs[0].j < n:
                sl2 = self._find_stem_matches(sl1)

                if sl2:
                    inequality = gp.LinExpr([2,-1,1],[sl1.base_pairs[0].var, sl1.var, sl2.var])
                    self.model.addConstr(inequality <= 1, f'FPIT-{sl1.base_pairs[0].i}-{sl1.base_pairs[0].j}')
                    inequality = gp.LinExpr([1,-1,-1],[sl1.var, sl2.var, sl1.base_pairs[0].var])
                    self.model.addConstr(inequality <= 0, f'FPOI-{sl1.base_pairs[0].i}-{sl1.base_pairs[0].j}')
            else:
                inequality = gp.LinExpr([2,-1],[sl1.base_pairs[0].var, sl1.var])
                self.model.addConstr(inequality <= 1, f'FPIT-{sl1.base_pairs[0].i}-{sl1.base_pairs[0].j}')
                inequality = gp.LinExpr([1,-1],[sl1.var, sl1.base_pairs[0].var])
                self.model.addConstr(inequality <= 0, f'FPOI-{sl1.base_pairs[0].i}-{sl1.base_pairs[0].j}')            
        self.model.update()

    def create_last_pair_constaints(self) -> None:
        n = len(self.rna_seq)

        for sl2 in self.stem_loops:
            if sl2.base_pairs[0].i > 1 and sl2.base_pairs[0].j < n:
                sl1 = self._find_stem_matches(sl2)

                if sl1:
                    inequality = gp.LinExpr([2,-1,1],[sl1.base_pairs[1].var, sl1.var, sl2.var])
                    self.model.addConstr(inequality <= 1, f'LPIT-{sl1.base_pairs[1].i}-{sl1.base_pairs[1].j}')
                    inequality = gp.LinExpr([1,-1,-1],[sl1.var, sl2.var, sl1.base_pairs[1].var])
                    self.model.addConstr(inequality <= 0, f'LPOI-{sl1.base_pairs[1].i}-{sl1.base_pairs[1].j}')
        self.model.update()

    def create_unpaired_nucleotides_constraints(self) -> None:
        n = len(self.rna_seq)

        for i in range(1, n + 1):
            inequality = gp.LinExpr(0)
            matches = self._find_base_pairs_with_index(self.base_pairs, i)

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


rna = seq_data['sequence']
print(rna)
rna_model = LILPModel(rna)
rna_model.create_base_pairs()
rna_model.create_hairpin_loops()
rna_model.create_stem_loops()
rna_model.create_internal_loops()
rna_model.create_bulge_loops()
rna_model.create_multi_loops()
rna_model.create_single_pair_constraints()
rna_model.create_no_crossing_constraints()
rna_model.create_stem_constraints()
rna_model.create_first_pair_constraints()
rna_model.create_last_pair_constaints()
rna_model.create_nucleotides()
rna_model.create_unpaired_nucleotides_constraints()
rna_model.create_hairpin_size_constraints()
rna_model.create_hairpin_ifthen_constraints()
rna_model.create_hairpin_onlyif_constraints()


# counter = 0

# print(type(rna_model.base_pairs))

# for bp in rna_model.base_pairs:    
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

dloop = RNALoop([rna_model.base_pairs[67], rna_model.base_pairs[90]],rna)

# P :: 113 :::: Q :: 24 :::: H :: 113 :::: I :: 1053 :::: B :: 387 :::: M :: 2694

