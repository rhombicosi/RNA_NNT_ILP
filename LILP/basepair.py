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
            
    def add_variable(self, model: gp.Model, label: str):
        if self.is_valid():
            name = f'{label}_{self.i}_{self.j}'
            self.var = model.addVar(vtype=GRB.BINARY, name=name)
        return self.var
 