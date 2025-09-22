import gurobipy as gp
from gurobipy import GRB
from utils.prepro_utils import *
from utils.sol_converter import *
from lilp_config import *
from basepair import *
from dloop import *
from stemloop import *
from hairpinloop import *
from internalloop import *
from bulgeloop import *
from multiloop import *

class LILP:
    def __init__(self, rna_seq: str, name: str):
        self.rna_seq = rna_seq
        self.model = gp.Model(name)
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

    def create_hairpin_vars(self) -> None:        
        # sorted_hairpins = sorted(self.hairpin_loops, key=lambda x: (x.size, x.energy))
        for hl in self.hairpin_loops:
            hl.add_variable(self.model)

    def create_stem_loops(self) -> None:
        for bp1 in self.base_pairs:
            bp2 = BasePair._find_base_pairs_matches(self.base_pairs, bp1.i + 1, bp1.j - 1)
            if bp2:
                stem = StemLoop([bp1, bp2], self.rna_seq)
                stem.add_variable(self.model)
                self.stem_loops.append(stem)
        self.model.update()

    def create_stem_vars(self) -> None:
        sorted_stems = sorted(self.stem_loops, key=lambda x: x.distance, reverse=True)   
        for sl in sorted_stems:
            sl.add_variable(self.model)

    def create_internal_loops(self) -> None:
        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                if bp2.i > bp1.i + 1 and bp2.j < bp1.j - 1:
                    internal = InternalLoop([bp1, bp2], self.rna_seq)
                    internal.add_variable(self.model)
                    self.internal_loops.append(internal)
        self.model.update()

    def create_internal_vars(self) -> None:
        sorted_internals = sorted(self.internal_loops, key=lambda x: (x.size, x.energy))
        for il in sorted_internals:
            il.add_variable(self.model)

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

    def add_stem_ifthen_constraints(self) -> None:
        for sl in self.stem_loops:
            sl.create_stem_ifthen_constraint(self.model)
        self.model.update()

    def add_stem_onlyif_constraints(self) -> None:
        for sl in self.stem_loops:
            sl.create_stem_onlyif_constraint(self.model)
        self.model.update()

    def add_stem_constraints(self) -> None:
        for sl in self.stem_loops:
            # if sl.energy > 0:
            #     sl.create_stem_ifthen_constraint(self.model)
            # else:
            #     sl.create_stem_onlyif_constraint(self.model)
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

    def add_hairpin_constraints(self) -> None:
        for hl in self.hairpin_loops: 
            if hl.energy > 0:
                hl.create_hairpin_ifthen_constraint(self.model, self.nucleotides)
                hl.create_hairpin_onlyif_constraint(self.model, self.base_pairs)
            else:
                hl.create_hairpin_onlyif_constraint(self.model, self.base_pairs)

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
    
    def add_internal_constraints(self)-> None:        
        for il in self.internal_loops:                
            if il.energy > 0:
                il.create_internal_ifthen_constraint(self.model, self.nucleotides)
                il.create_internal_onlyif_constraint(self.model, self.base_pairs)
            else:
                il.create_internal_onlyif_constraint(self.model, self.base_pairs)

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

    def add_bulge_constraints(self) -> None:
        for bl in self.bulge_loops:
            if bl.energy > 0:
                bl.create_bulge_ifthen_constraint(self.model, self.nucleotides)
            else:
                bl.create_bulge_onlyif_constraint(self.model, self.base_pairs)

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
    
    def add_multi_constraints(self) -> None:
        for ml in self.multi_loops:
            if ml.energy > 0:
                ml.create_multi_ifthen_constraint(self.model, self.nucleotides)
            else:
                ml.create_multi_onlyif_constraint(self.model, self.base_pairs)
        self.model.update()

    def add_multi_max_number_constraint(self) -> None:
        MultiLoop.create_multi_max_number_constraint(self.model, self.multi_loops)

    def create_stem_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([sl.energy for sl in self.stem_loops], [sl.var for sl in self.stem_loops])
        return objective
    
    def create_hairpin_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([hl.energy for hl in self.hairpin_loops], [hl.var for hl in self.hairpin_loops])
        return objective
    
    def create_internal_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([il.energy for il in self.internal_loops], [il.var for il in self.internal_loops])
        return objective
    
    def create_bulge_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([bl.energy for bl in self.bulge_loops], [bl.var for bl in self.bulge_loops])
        return objective
    
    def create_multi_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([ml.energy for ml in self.multi_loops], [ml.var for ml in self.multi_loops])
        return objective
    
    def create_objective(self, stem, hairpin, internal, bulge, multi) -> gp.LinExpr:
        objective = gp.LinExpr()
        if stem:
            objective.add(self.create_stem_term())
        if hairpin:
            objective.add(self.create_hairpin_term())
        if internal:
            objective.add(self.create_internal_term())
        if bulge:
            objective.add(self.create_bulge_term())
        if multi:
            objective.add(self.create_multi_term())
        self.model.setObjective(objective, GRB.MINIMIZE)

    def create_variables(self, stem, hairpin, internal, bulge, multi):
        self.create_base_pairs()       
        if hairpin:
            self.create_hairpin_loops() 
        if stem:            
            self.create_stem_loops()                                
        self.create_nucleotides()
        if internal:  
            self.create_internal_loops()  
        if bulge:
            self.create_bulge_loops()
        if multi:
            self.create_multi_loops()
    
    def create_constraints(self, stem, hairpin, internal, bulge, multi):
        self.add_single_pair_constraints()
        self.add_no_crossing_constraints()
        if stem:
            # self.add_stem_ifthen_constraints()
            # self.add_stem_onlyif_constraints()
            self.add_stem_constraints()
        self.add_unpaired_nucleotides_constraints()
        if hairpin:
            self.add_hairpin_size_constraints()
            self.add_hairpin_constraints()
            # self.add_hairpin_ifthen_constraints()
            # self.add_hairpin_onlyif_constraints()            
            self.add_hairpin_max_number_constraint()
        if internal:
            self.add_internal_size_constraints()
            self.add_internal_constraints()
            # self.add_internal_onlyif_constraints()
            # self.add_internal_ifthen_constraints()            
            self.add_internal_max_number_constraint()
            # self.model.addConstr(self.model.getVarByName(f'INTERNAL_5_39_11_36') == 1)
        if bulge:
            self.add_bulge_size_constraints()
            self.add_bulge_constraints()
            # self.add_bulge_ifthen_constraints()
            # self.add_bulge_onlyif_constraints()
            # self.add_bulge_max_number_constraint()
        if multi:
            self.add_multi_size_constraints()
            self.add_multi_constraints()
            # self.add_multi_ifthen_constraints()
            # self.add_multi_onlyif_constraints()
            # self.add_multi_max_number_constraint()
        