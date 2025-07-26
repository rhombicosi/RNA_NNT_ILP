from dloop import *

class HairpinLoop(Loop):

    def __init__(self, base_pairs, RNA):
        super().__init__(base_pairs, RNA)

    def calculate_energy() -> int:
        pass

    def create_hairpin_size_constraint(self, model: gp.Model) -> None:
        if not self.is_valid_size():
            inequality = gp.LinExpr([1], [self.var])
            model.addConstr(inequality == 0, f'HS-{self.base_pairs[0].i}-{self.base_pairs[0].j}')

    def create_hairpin_ifthen_constraint(self, model: gp.Model, nucleotides: List[gp.Var]) -> None:
        inequality = gp.LinExpr(0)
        for u in range(self.base_pairs[0].i + 1, self.base_pairs[0].j):
            inequality.add(gp.LinExpr([1], [nucleotides[u - 1]]))
        
        inequality.add(gp.LinExpr([1, -1],[self.base_pairs[0].var, self.var]))            
        model.addConstr(inequality <= self.size, f'HIT-{self.base_pairs[0].i}-{self.base_pairs[0].j}')

    def create_hairpin_onlyif_constraint(self, model: gp.Model, base_pairs: List[BasePair]) -> None:        
        for u in range(self.base_pairs[0].i + 1, self.base_pairs[0].j):
            inequality = gp.LinExpr([2], [self.var])
            matches = BasePair._find_base_pairs_with_index(base_pairs, u)

            for bp in matches:
                inequality.add(gp.LinExpr([1], [bp.var]))
            
            inequality.add(gp.LinExpr([-1], [self.base_pairs[0].var]))
            model.addConstr(inequality <= 1, f'HOI-{self.base_pairs[0].i}-{self.base_pairs[0].j}-{u}')

    def create_hairpin_max_number_constraint(model: gp.Model, hairpin_loops: List["HairpinLoop"]) -> None:
        inequality = gp.LinExpr(0)

        for hl in hairpin_loops:
            inequality.add(gp.LinExpr([1], [hl.var]))
        model.addConstr(inequality <= MAX_NUM_OF_LOOPS[hl.type], f'HMN')
        model.update()    
