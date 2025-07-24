from typing import List
import gurobipy as gp
from gurobipy import GRB
from lilp_config import *

class BasePair:

    def __init__(self, i: int, j: int, rna: str):
        self.i = i
        self.j = j
        self.rna = rna
        self.nt1 = rna[i-1]
        self.nt2 = rna[j-1]
        self.var = None

    def distance(self) -> int:
        return self.j - self.i

    def is_valid(self) -> bool:
        return f'{self.nt1}{self.nt2}' in VALID_PAIRS and BasePair.distance(self) > MIN_D
            
    def add_variable(self, model: gp.Model, label: str) -> gp.Var:
        if self.is_valid():
            name = f'{label}_{self.i}_{self.j}'
            self.var = model.addVar(vtype=GRB.BINARY, name=name)
        return self.var
    
    def _find_base_pairs_with_index(base_pairs: List["BasePair"], index: int) -> List["BasePair"]:
        return [bp for bp in base_pairs if bp.i == index or bp.j == index]
    
    def _find_base_pairs_matches(base_pairs: List["BasePair"], i: int, j: int) -> "BasePair":
        return next((bp for bp in base_pairs if bp.i == i and bp.j == j), None)
    
    def create_single_pair_constraint(model: gp.Model, base_pairs: List["BasePair"], i: int) -> None:        
        inequality = gp.LinExpr(0)
        matches = BasePair._find_base_pairs_with_index(base_pairs, i)

        if matches:
            for bp in matches:
                inequality.add(gp.LinExpr([1.0], [bp.var]))
            model.addConstr(inequality <= 1, f'SP-{i}')

    def create_no_crossing_constraint(model: gp.Model, bp1: "BasePair", bp2: "BasePair") -> None:        
        if bp2.i > bp1.i and bp2.i < bp1.j and bp2.j > bp1.j:
            inequality = gp.LinExpr(0)
            inequality.add(gp.LinExpr([1.0, 1.0],[bp1.var, bp2. var]))
            model.addConstr(inequality <= 1, f'NC-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}')
