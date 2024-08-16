from binary_variables_grb import *
from constraints_grb import *
from objective_grb import *
from pathlib import Path
import time
import prepro_run

try:
    mip.setObjective(objectiveTerm(RNA), GRB.MINIMIZE)
    # mip.addConstr(objectiveTerm(RNA) >= MFE,"CMFE")
    # mip.addConstr(mip.getVarByName(f'B(14,15,28,30)') == 1)
    # mip.addConstr(mip.getVarByName(f'H(19,24)') == 1)
    onePairConstraints(RNA)
    noCrossConstraints(RNA)
    stemConstraints(RNA)
    firstPairConstraints(RNA)
    lastPairConstraints(RNA)
    hairpinNTConstraints(RNA)
    hairpinIfThenConstraints(RNA)
    hairpinOnlyIfConstraints(RNA)
    numHairpinConstraints(RNA, numH)
    internalNTConstraints(RNA)
    internalIfThenConstraints(RNA)
    internalOnlyIfConstraints(RNA)
    numInternalConstraints(RNA, numI)
    bulgeNTConstraints(RNA)
    bulgeIfThenConstraints(RNA)
    bulgeOnlyIfConstraints(RNA)
    numBulgeConstraints(RNA, numB)

    mip.update()
    
    lp_file_name = seq_data["identifier"].split(" ")[0]

    start_time = time.time()
    mip.write(f'{prepro_run.lp_dir}/{lp_file_name}-decomposition-grb.lp')
    print(f"--- {time.time() - start_time} seconds ---")

    # uncomment optimization tuning params if needed
    # mip.setParam("MIPGap", 0.2)
    # mip.setParam("PoolSolutions", 20)
    # mip.setParam("PoolSearchMode", 2)
    # mip.setParam("SolFiles", f"{chain_f}-decomposition-grb")

    mip.optimize()

    mip.write(f'{prepro_run.sol_dir}/{lp_file_name}-decomposition-grb.sol')

    print(f'Obj: {mip.ObjVal:g}')

except gp.GurobiError as e:
    print(f'Error code {e.errno}: {e}')

except AttributeError:
    print('Encountered an attribute error')
