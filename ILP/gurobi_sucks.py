import gurobipy as gp


m = gp.read("3G4S_9-decomposition.lp")
m.optimize()
m.write("3G4S_9.sol")