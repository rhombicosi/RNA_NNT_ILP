from dloop import *
from lilp_config import *

class MultiLoop(Loop):

    def __init__(self, base_pairs, RNA):
        super().__init__(base_pairs, RNA)
        self.energy = self.calculate_energy()

    def calculate_energy(self) -> int:
        if self.is_valid_size():
            G = SCALE * (B * self.size + C * self.degree)
        else:
            G = M
        return round(G)

    def create_multi_size_constraint(self, model: gp.Model) -> None: 
        bp1 = self.base_pairs[0]
        bp2 = self.base_pairs[1]
        bp3 = self.base_pairs[2]

        if not self.is_valid_size():
            inequality = gp.LinExpr([1], [self.var])
            model.addConstr(inequality == 0, f'MS-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}-{bp3.i}-{bp3.j}')
    
    def create_multi_ifthen_constraint(self, model: gp.Model, nucleotides: List[gp.Var]) -> None:
        bp1 = self.base_pairs[0]
        bp2 = self.base_pairs[1]
        bp3 = self.base_pairs[2]

        inequality = gp.LinExpr(0)                       

        for u in range(bp1.i + 1, bp2.i):
            inequality.add(gp.LinExpr([1],[nucleotides[u - 1]]))

        for u in range(bp2.j + 1, bp3.i):
            inequality.add(gp.LinExpr([1],[nucleotides[u - 1]]))
        
        for u in range(bp3.j + 1, bp1.j):
            inequality.add(gp.LinExpr([1],[nucleotides[u - 1]]))
            
        inequality.add(gp.LinExpr([1, 1, 1, -1],[bp1.var, bp2.var, bp3.var, self.var]))
        model.addConstr(inequality <= self.size + 1, f'MIT-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}-{bp3.i}-{bp3.j}')

    def create_multi_onlyif_constraint(self, model: gp.Model, base_pairs: List[BasePair]) -> None:
        bp1 = self.base_pairs[0]
        bp2 = self.base_pairs[1]
        bp3 = self.base_pairs[2]

        if bp1.i + 1 == bp2.i or bp2.j + 1 == bp3.i or bp3.j + 1 == bp1.j:
            inequality = gp.LinExpr(0)
            inequality.add(gp.LinExpr([4], [self.var]))
            inequality.add(gp.LinExpr([-1, -1, -1], [bp1.var, bp2.var, bp3.var]))
            model.addConstr(inequality <= 1, f'MOI-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}-{bp3.i}-{bp3.j}')
        elif bp1.i + 1 != bp2.i:
            for u in range(bp1.i + 1, bp2.i):            
                inequality = gp.LinExpr(0)
                inequality.add(gp.LinExpr([4], [self.var]))

                matches = BasePair._find_base_pairs_with_index(base_pairs, u)
                
                for bp in matches:
                    inequality.add(gp.LinExpr([1], [bp.var]))
                
                inequality.add(gp.LinExpr([-1, -1, -1],[bp1.var, bp2.var, bp3.var]))
                model.addConstr(inequality <= 1, f'MOI-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}-{bp3.i}-{bp3.j}-{u}')        
        elif bp2.j + 1 != bp3.i:
            for u in range(bp2.j + 1, bp3.i):
                inequality = gp.LinExpr(0)
                inequality.add(gp.LinExpr([4], [self.var]))

                matches = BasePair._find_base_pairs_with_index(base_pairs, u)
                
                for bp in matches:
                    inequality.add(gp.LinExpr([1], [bp.var]))
                
                inequality.add(gp.LinExpr([-1, -1, -1],[bp1.var, bp2.var, bp3.var]))
                model.addConstr(inequality <= 1, f'MOI-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}-{bp3.i}-{bp3.j}-{u}') 
        elif bp3.j + 1 != bp1.j:
            for u in range(bp3.j + 1, bp1.j):
                inequality = gp.LinExpr(0)
                inequality.add(gp.LinExpr([4], [self.var]))

                matches = BasePair._find_base_pairs_with_index(base_pairs, u)
                
                for bp in matches:
                    inequality.add(gp.LinExpr([1], [bp.var]))
                
                inequality.add(gp.LinExpr([-1, -1, -1],[bp1.var, bp2.var, bp3.var]))
                model.addConstr(inequality <= 1, f'MOI-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}-{bp3.i}-{bp3.j}-{u}')

    def create_multi_max_number_constraint(model: gp.Model, multi_loops: List["MultiLoop"]) -> None:
        inequality = gp.LinExpr(0)

        for ml in multi_loops:
            inequality.add(gp.LinExpr([1], [ml.var]))
        model.addConstr(inequality <= MAX_NUM_OF_LOOPS[ml.type], f'MMN')
        model.update 

# RNA = 'GCCGCGAACCCCGCCAGGCCCGGAAGGGAGCAACGGUAGUGGUGGAU'
# bp1 = BasePair(10,35,RNA)
# bp2 = BasePair(11,18,RNA)
# bp3 = BasePair(19,28,RNA)
# multi = MultiLoop((bp1,bp2,bp3), RNA)
# print(multi.size)
# print(multi.energy)

