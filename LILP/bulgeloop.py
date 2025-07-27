from dloop import *

class BulgeLoop(Loop):

    def __init__(self, base_pairs, RNA):
        super().__init__(base_pairs, RNA)

    def calculate_energy() -> int:
        pass

    def create_bulge_size_constraint(self, model: gp.Model) -> None:        
        if not self.is_valid_size():
            inequality = gp.LinExpr([1], [self.var])
            model.addConstr(inequality == 0, f'BS-{self.base_pairs[0].i}-{self.base_pairs[0].j}-{self.base_pairs[1].i}-{self.base_pairs[1].j}')

    def create_bulge_ifthen_constraint(self, model: gp.Model, nucleotides: List[gp.Var]) -> None:
        bp1 = self.base_pairs[0]
        bp2 = self.base_pairs[1]  

        if bp2.i == bp1.i + 1:
            inequality = gp.LinExpr(0)   
        
            for u in range(bp2.j + 1, bp1.j):
                inequality.add(gp.LinExpr([1], [nucleotides[u - 1]]))
            inequality.add(gp.LinExpr([1, 1, -1], [bp1.var, bp2.var, self.var]))
            model.addConstr(inequality <= self.size + 1, f'BIT-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}')
        elif bp2.j == bp1.j - 1:
            inequality = gp.LinExpr(0)

            for u in range(bp1.i + 1, bp2.i):
                inequality.add(gp.LinExpr([1], [nucleotides[u - 1]]))
            inequality.add(gp.LinExpr([1, 1, -1], [bp1.var, bp2.var, self.var]))
            model.addConstr(inequality <= self.size + 1, f'BIT-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}')

    def create_bulge_onlyif_constraint(self, model: gp.Model, base_pairs: List[BasePair]) -> None:
        bp1 = self.base_pairs[0]
        bp2 = self.base_pairs[1] 

        if bp2.i == bp1.i + 1:
            for u in range(bp2.j + 1, bp1.j):
                inequality = gp.LinExpr(0)
                inequality.add(gp.LinExpr([3], [self.var]))

                matches = BasePair._find_base_pairs_with_index(base_pairs, u)
                
                for bp in matches:
                    inequality.add(gp.LinExpr([1], [bp.var]))

                inequality.add(gp.LinExpr([-1, -1],[bp1.var, bp2.var]))
                model.addConstr(inequality <= 1, f'BOI-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}-{u}')        
        elif bp2.j == bp1.j - 1:
            inequality = gp.LinExpr(0)

            for u in range(bp1.i + 1, bp2.i):
                inequality = gp.LinExpr(0)
                inequality.add(gp.LinExpr([3], [self.var]))

                matches = BasePair._find_base_pairs_with_index(base_pairs, u)
                
                for bp in matches:
                    inequality.add(gp.LinExpr([1], [bp.var]))

                inequality.add(gp.LinExpr([-1, -1],[bp1.var, bp2.var]))
                model.addConstr(inequality <= 1, f'BOI-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}-{u}') 

    def create_bulge_max_number_constraint(model: gp.Model, bulge_loops: List["BulgeLoop"]) -> None:
        inequality = gp.LinExpr(0)

        for bl in bulge_loops:
            inequality.add(gp.LinExpr([1], [bl.var]))
        model.addConstr(inequality <= MAX_NUM_OF_LOOPS[bl.type], f'BMN')
        model.update
