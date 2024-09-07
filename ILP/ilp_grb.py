from binary_variables_grb import *
from constraints_grb import *
from objective_grb import *
from pathlib import Path
import time
from prepro_run import *


def optimize(numH, numI, numB, seq_files, seq_number, lp_dir, sol_dir):

    
    chain_file = seq_files[seq_number]
    chain_name_with_ext = os.path.basename(chain_file)        
    chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
    lp_file_name = chain_name_without_ext
    seq_data = parse_seq_file(chain_file)
    this_RNA = seq_data['sequence']

    mip = gp.Model(f'MIP-{seq_number}')

    try:
        
        listP, listQ, listF, listL, listH, listI, listB, listX, listY, listZ = add_binary_vars(this_RNA, mip)

        mip.setObjective(objectiveTerm(this_RNA, listQ, listH, listI, listB), GRB.MINIMIZE)

        
        # mip.addConstr(objectiveTerm(RNA) >= MFE,"CMFE")
        # mip.addConstr(mip.getVarByName(f'B(14,15,28,30)') == 1)
        # mip.addConstr(mip.getVarByName(f'H(19,24)') == 1)
        onePairConstraints(this_RNA, mip)
        noCrossConstraints(this_RNA, mip)
        stemConstraints(this_RNA, mip)
        firstPairConstraints(this_RNA, mip)
        lastPairConstraints(this_RNA, mip)
        hairpinNTConstraints(this_RNA, mip)
        hairpinIfThenConstraints(this_RNA, mip)
        hairpinOnlyIfConstraints(this_RNA, mip)
        numHairpinConstraints(this_RNA, numH, mip)
        internalNTConstraints(this_RNA, mip)
        internalIfThenConstraints(this_RNA, mip)
        internalOnlyIfConstraints(this_RNA, mip)
        numInternalConstraints(this_RNA, numI, mip)
        bulgeNTConstraints(this_RNA, mip)
        bulgeIfThenConstraints(this_RNA, mip)
        bulgeOnlyIfConstraints(this_RNA, mip)
        numBulgeConstraints(this_RNA, numB, mip)

        mip.update()   

        start_time = time.time()
        mip.write(f'{lp_dir}/{lp_file_name}-loopdeco.lp')
        print(f"--- {time.time() - start_time} seconds ---")

        # uncomment optimization tuning params if needed
        # mip.setParam("MIPGap", 0.2)
        # mip.setParam("PoolSolutions", 20)
        # mip.setParam("PoolSearchMode", 2)
        # mip.setParam("SolFiles", f"{chain_f}-decomposition-grb")

        mip.optimize()

        mip.write(f'{sol_dir}/{lp_file_name}-loopdeco.sol')

        print(f'Obj: {mip.ObjVal:g}')
        
        obj_val = mip.ObjVal

    except gp.GurobiError as e:
        print(f'Error code {e.errno}: {e}')
        obj_val = None

    except AttributeError:
        print('Encountered an attribute error')
        obj_val = None

    if obj_val is not None:
        return obj_val, lp_file_name
    else:
        print("Object value was not assigned due to an error.")
