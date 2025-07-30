import os
from pathlib import Path
from utils.prepro_utils import *
import gurobipy as gp
from gurobipy import GRB
from lilp_config import *
from basepair import *
from dloop import *
from stemloop import *
from hairpinloop import *
from internalloop import *
from bulgeloop import *
from multiloop import *

class LILPModel:
    def __init__(self, rna_seq: str):
        self.rna_seq = rna_seq
        self.model = gp.Model(f'MIP')
        self.nucleotides : List[gp.Var] = []
        self.base_pairs : List[BasePair] = []
        self.first_pairs : List[BasePair] = []
        self.last_pairs : List[BasePair] = []
        self.hairpin_loops : List[HairpinLoop] = []
        self.stem_loops : List[StemLoop] = []
        self.internal_loops : List[InternalLoop] = []
        self.bulge_loops : List[BulgeLoop] = []
        self.multi_loops : List[MultiLoop] = []
        self.objective : gp.LinExpr   

    def create_nucleotides(self) -> None:
        n = len(self.rna_seq)

        for i in range(1, n + 1):
            var = self.model.addVar(vtype=GRB.BINARY, name=f'X_{i}')   
            self.nucleotides.append(var)     
        self.model.update()

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

    def create_hairpin_loops(self) -> None:

        for bp in self.base_pairs:
            hairpin = HairpinLoop([bp], self.rna_seq)
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
                    internal = InternalLoop([bp1, bp2], self.rna_seq)
                    internal.add_variable(self.model)
                    self.internal_loops.append(internal)
        self.model.update()

    def create_bulge_loops(self) -> None:

        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                if (bp2.i == bp1.i + 1 and bp2.j < bp1.j - 1) or (bp2.i > bp1.i + 1 and bp2.j == bp1.j - 1):
                    bulge = BulgeLoop([bp1, bp2], self.rna_seq)
                    bulge.add_variable(self.model)
                    self.bulge_loops.append(bulge)
        self.model.update()

    def create_multi_loops(self) -> None:

        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                for bp3 in self.base_pairs:
                    if bp2.i > bp1.i and bp3.i > bp2.j and bp1.j > bp3.j:
                        multi = MultiLoop([bp1, bp2, bp3], self.rna_seq)
                        multi.add_variable(self.model)
                        self.multi_loops.append(multi)
        self.model.update()

    def add_unpaired_nucleotides_constraints(self) -> None:
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

    def add_hairpin_size_constraints(self) -> None:
        for hl in self.hairpin_loops:
            hl.create_hairpin_size_constraint(self.model)
        self.model.update()

    def add_hairpin_ifthen_constraints(self) -> None:
        for hl in self.hairpin_loops:
            hl.create_hairpin_ifthen_constraint(self.model, self.nucleotides)
        self.model.update()

    def add_hairpin_onlyif_constraints(self) -> None:
        for hl in self.hairpin_loops:
            hl.create_hairpin_onlyif_constraint(self.model, self.base_pairs)
        self.model.update()

    def add_hairpin_max_number_constraint(self) -> None:
        HairpinLoop.create_hairpin_max_number_constraint(self.model, self.hairpin_loops)

    def add_internal_size_constraints(self) -> None:
        for il in self.internal_loops:
            il.create_internal_size_constraint(self.model)
        self.model.update()

    def add_internal_ifthen_constraints(self) -> None:
        for il in self.internal_loops:
            il.create_internal_ifthen_constraint(self.model, self.nucleotides)
        self.model.update()

    def add_internal_onlyif_constraints(self) -> None:
        for il in self.internal_loops:
            il.create_internal_onlyif_constraint(self.model, self.base_pairs)
        self.model.update()

    def add_internal_max_number_constraint(self) -> None:
        InternalLoop.create_internal_max_number_constraint(self.model, self.internal_loops)

    def add_bulge_size_constraints(self) -> None:
        for il in self.bulge_loops:
            il.create_bulge_size_constraint(self.model)
        self.model.update()

    def add_bulge_ifthen_constraints(self) -> None:
        for bl in self.bulge_loops:
            bl.create_bulge_ifthen_constraint(self.model, self.nucleotides)
        self.model.update()

    def add_bulge_onlyif_constraints(self) -> None:
        for bl in self.bulge_loops:
            bl.create_bulge_onlyif_constraint(self.model, self.base_pairs)
        self.model.update()

    def add_bulge_max_number_constraint(self) -> None:
        BulgeLoop.create_bulge_max_number_constraint(self.model, self.bulge_loops)

    def add_multi_size_constraints(self) -> None:
        for ml in self.multi_loops:
            ml.create_multi_size_constraint(self.model)
        self.model.update()

    def add_multi_ifthen_constraints(self) -> None:
        for ml in self.multi_loops:
            ml.create_multi_ifthen_constraint(self.model, self.nucleotides)
        self.model.update()

    def add_multi_onlyif_constraints(self) -> None:
        for ml in self.multi_loops:
            ml.create_multi_onlyif_constraint(self.model, self.base_pairs)
        self.model.update()

    def add_multi_max_number_constraint(self) -> None:
        MultiLoop.create_multi_max_number_constraint(self.model, self.multi_loops)

    def create_stem_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([sl.energy for sl in self.stem_loops], [sl.var for sl in self.stem_loops])
        return objective
    
    def create_first_pair_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([fp.pair_penalty_energy for fp in self.first_pairs], [fp.var for fp in self.first_pairs])
        return objective


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
rna_model.create_internal_loops()
rna_model.create_bulge_loops()
rna_model.create_multi_loops()

rna_model.add_single_pair_constraints()
rna_model.add_no_crossing_constraints()

rna_model.add_stem_constraints()
rna_model.add_first_pair_constraints()
rna_model.add_last_pair_constraints()

rna_model.create_nucleotides()
rna_model.add_unpaired_nucleotides_constraints()

rna_model.add_hairpin_size_constraints()
rna_model.add_hairpin_ifthen_constraints()
rna_model.add_hairpin_onlyif_constraints()
rna_model.add_hairpin_max_number_constraint()

rna_model.add_internal_size_constraints()
rna_model.add_internal_ifthen_constraints()
rna_model.add_internal_onlyif_constraints()
rna_model.add_internal_max_number_constraint()

rna_model.add_bulge_size_constraints()
rna_model.add_bulge_ifthen_constraints()
rna_model.add_bulge_onlyif_constraints()
rna_model.add_bulge_max_number_constraint()

rna_model.add_multi_size_constraints()
rna_model.add_multi_ifthen_constraints()
rna_model.add_multi_onlyif_constraints()
rna_model.add_multi_max_number_constraint()

rna_model.model.setObjective(rna_model.create_stem_term(), GRB.MINIMIZE)

rna_model.model.write(f'lilp-{seq_no}.lp')

counter = 0

# print(type(rna_model.base_pairs))

# for bp in rna_model.base_pairs:    
#     print(f"[{counter}]::Base pair: ({bp.i}, {bp.j}), First pair {bp.nt1 + bp.nt2} energy: {bp.pair_penalty_energy}, Variable name: {bp.var.VarName}, Gurobi var: {bp.var}")
#     counter+=1

# for bp in rna_model.first_pairs:    
#     print(f"[{counter}]::Base pair: ({bp.i}, {bp.j}), Variable name: {bp.var.VarName}, Gurobi var: {bp.var}")
#     counter+=1

for h in rna_model.hairpin_loops:    
    print(f"[{counter}]:: Variable name: {h.var.VarName}, Energy: {h.energy}, Gurobi var: {h.var}")
    counter+=1

# for s in rna_model.stem_loops:    
#     print(f"[{counter}]:: Variable name: {s.var.VarName}, Energy of {rna[s.first_pair.i - 1] + rna[s.first_pair.j - 1]} followed by {rna[s.last_pair.i - 1] + rna[s.last_pair.j - 1]}: {s.energy}, Gurobi var: {s.var}")
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

