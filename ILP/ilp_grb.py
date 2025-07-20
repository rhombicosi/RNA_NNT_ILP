from binary_variables_grb import *
from constraints_grb import *
from objective_grb import *
from pathlib import Path
import time
from prepro_run import *


def optimize(seq_files, seq_number, lp_dir, sol_dir, s_start, start = solstart_dir):
    
    chain_file = seq_files[seq_number]
    chain_name_with_ext = os.path.basename(chain_file)        
    chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
    lp_file_name = chain_name_without_ext
    seq_data = parse_seq_file(chain_file)
    this_RNA = seq_data['sequence']

    print(chain_name_without_ext)
    print(this_RNA)

    mip = gp.Model(f'MIP-{seq_number}')

    # set the number of loops
    numH = 1 #len(this_RNA)//5
    numI = 1 #len(this_RNA)//4
    numB = 0 #len(this_RNA)//3
    numM = 1

    try:
        
        listP, listQ, listH, listI, listB, listM, listX = add_binary_vars(this_RNA, mip, 0)

        mip.setObjective(objectiveTerm(this_RNA, listQ, listH, listI, listB, listM), GRB.MINIMIZE)
        
        # mip.addConstr(objectiveTerm(this_RNA, listQ, listH, listI, listB, listM) >= MFE,"CMFE")
        # mip.addConstr(mip.getVarByName(f'B(14,15,28,30)') == 1)
        onePairConstraints(this_RNA, mip)
        noCrossConstraints(this_RNA, mip)
        stemConstraints(this_RNA, mip)
        # firstPairConstraints(this_RNA, mip)
        # lastPairConstraints(this_RNA, mip)
        # hairpinNTConstraints(this_RNA, mip)
        hairpinZeroConstraints(this_RNA, mip, maxH)
        hairpinIfThenConstraints(this_RNA, mip)
        hairpinOnlyIfConstraints(this_RNA, mip)
        # numHairpinConstraints(this_RNA, numH, mip)
        # internalNTConstraints(this_RNA, mip)
        internalZeroConstraints(this_RNA, mip, maxI)
        internalIfThenConstraints(this_RNA, mip)
        internalOnlyIfConstraints(this_RNA, mip)
        # numInternalConstraints(this_RNA, numI, mip)
        # bulgeNTConstraints(this_RNA, mip)
        bulgeZeroConstraints(this_RNA, mip, maxB)
        bulgeIfThenConstraints(this_RNA, mip)
        bulgeOnlyIfConstraints(this_RNA, mip)
        # numBulgeConstraints(this_RNA, numB, mip)
        # multiNTConstraints(this_RNA, mip)
        multiZeroConstraints(this_RNA,mip,maxM)
        multiIfThenConstraints(this_RNA, mip)
        multiOnlyIfConstraints(this_RNA, mip)
        numMultiConstraints(this_RNA, numM, mip)
        loopNTConstraints(this_RNA, mip)

        mip.update()   

        start_time = time.time()
        mip.write(f'{lp_dir}/{lp_file_name}-loopdeco.lp')
        print(f"--- {time.time() - start_time} seconds ---")

        # uncomment optimization tuning params if needed
        # mip.setParam("MIPGap", 0.2)
        # mip.setParam("PoolSolutions", 20)
        # mip.setParam("PoolSearchMode", 2)
        # mip.setParam("SolFiles", f"{chain_f}-decomposition-grb")
        mip.setParam("LogFile", f'{grb_log_dir}/log-{lp_file_name}')
        # mip.setParam("Presolve", 2)
        # mip.setParam("Cuts", 3)
        # mip.setParam("CliqueCuts", 2)
        # mip.setParam("MIPFocus", 1)

        if start:

            mip.NumStart = s_start
            mip.update()

            # iterate over all MIP starts
            for s in range(0,mip.NumStart):
            
                # set StartNumber
                mip.params.StartNumber = s

                solvars = read_sol(f'{solstart_dir}/{lp_file_name}-start-{s}.sol')

                # now set MIP start values using the Start attribute, e.g.:
                for v in mip.getVars(): 
                    if v.varName in solvars.keys():         
                        v.Start = round(int(solvars[v.varName]), 1)            
                    
            mip.update()

        vars_with_Q = [var for var in mip.getVars() if 'Q' in var.VarName]
        for i in range(len(vars_with_Q)):
            vars_with_Q[i].setAttr("BranchPriority",100)

        # vars_with_H = [var for var in mip.getVars() if 'H' in var.VarName]
        # for i in range(len(vars_with_H)):
        #     vars_with_H[i].setAttr("BranchPriority",10)

        # vars_with_I = [var for var in mip.getVars() if 'I' in var.VarName]
        # for i in range(len(vars_with_I)):
        #     vars_with_I[i].setAttr("BranchPriority",100)

        # vars_with_B = [var for var in mip.getVars() if 'B' in var.VarName]
        # for i in range(len(vars_with_B)):
        #     vars_with_B[i].setAttr("BranchPriority",50)

        # vars_with_M = [var for var in mip.getVars() if 'M' in var.VarName]
        # for i in range(len(vars_with_M)):
        #     vars_with_M[i].setAttr("BranchPriority",50)

        mip.update()

        opt_start_time = time.time()
        mip.optimize()
        opt_time = time.time() - opt_start_time

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
        return obj_val, lp_file_name, opt_time
    else:
        print("Object value was not assigned due to an error.")

def optimize_start(seq_files, seq_number, lpstart_dir, solstart_dir):
    
    chain_file = seq_files[seq_number]
    chain_name_with_ext = os.path.basename(chain_file)        
    chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
    lp_file_name = chain_name_without_ext
    seq_data = parse_seq_file(chain_file)
    this_RNA = seq_data['sequence']

    print(chain_name_without_ext)
    print(this_RNA)

    mip = gp.Model(f'MIP-{seq_number}')

    # set the number of loops
    numH = 1#len(this_RNA)//5
    numI = 3#len(this_RNA)//4
    numB = 0#len(this_RNA)//3
    numM = 0#len(this_RNA)//3

    start = 1

    try:
        
        listP, listQ, listH, listI, listX = add_binary_vars(this_RNA, mip, start)
        mip.update()

        mip.setObjective(objectiveStartTerm(this_RNA, listQ, listH, listI), GRB.MINIMIZE)
        
        # mip.addConstr(objectiveStartTerm(this_RNA, listQ, listH, listI) >= MFE,"CMFE")
        # mip.addConstr(mip.getVarByName(f'B(14,15,28,30)') == 1)
        # mip.addConstr(mip.getVarByName(f'H(19,24)') == 1)
        # mip.addConstr(mip.getVarByName(f'I(6,12,37,41)') == 1)
        onePairConstraints(this_RNA, mip)
        noCrossConstraints(this_RNA, mip)
        stemConstraints(this_RNA, mip)
        # firstPairConstraints(this_RNA, mip)
        # lastPairConstraints(this_RNA, mip)
        # hairpinNTConstraints(this_RNA, mip)
        hairpinZeroConstraints(this_RNA, mip, maxH)
        hairpinIfThenConstraints(this_RNA, mip)
        # hairpinOnlyIfConstraints(this_RNA, mip)
        # numHairpinConstraints(this_RNA, numH, mip)
        # internalNTConstraints(this_RNA, mip)
        internalZeroConstraints(this_RNA, mip, maxI)
        internalIfThenConstraints(this_RNA, mip)
        internalOnlyIfConstraints(this_RNA, mip)
        # numInternalConstraints(this_RNA, numI, mip)
        # bulgeNTConstraints(this_RNA, mip)
        # bulgeZeroConstraints(this_RNA, mip, maxB)
        # bulgeIfThenConstraints(this_RNA, mip)
        # bulgeOnlyIfConstraints(this_RNA, mip)
        # numBulgeConstraints(this_RNA, numB, mip)
        # multiNTConstraints(this_RNA, mip)
        # multiZeroConstraints(this_RNA,mip,maxM)
        # multiIfThenConstraints(this_RNA, mip)
        # multiOnlyIfConstraints(this_RNA, mip)
        # numMultiConstraints(this_RNA, numM, mip)
        # numGreaterMultiConstraints(this_RNA, numM, mip)
        # consecutiveUnpairedConstraints(this_RNA, L, mip)
        loopNTConstraints(this_RNA, mip)

        mip.update()   

        start_time = time.time()
        mip.write(f'{lpstart_dir}/{lp_file_name}-start.lp')
        print(f"--- {time.time() - start_time} seconds ---")

        # uncomment optimization tuning params if needed
        # mip.setParam("MIPGap", 0.2)
        # mip.setParam("PoolSolutions", 20)
        # mip.setParam("PoolSearchMode", 2)
        # mip.setParam("SolFiles", f"{chain_f}-decomposition-grb")
        mip.setParam("LogFile", f'{grb_log_dir}/log-{lp_file_name}')
        # mip.setParam("MIPFocus", 2)
        # mip.setParam("ConcurrentMIP ", 5)
        # mip.setParam("Presolve", 2)
        # mip.setParam("Cuts", 3)
        # mip.setParam("CliqueCuts", 2)

        # vars_with_Q = [var for var in mip.getVars() if 'Q' in var.VarName]
        # for i in range(len(vars_with_Q)):
        #     vars_with_Q[i].setAttr("BranchPriority",100)
        
        
        # mip.update()
        
        opt_start_time = time.time()        
        mip.optimize()
        opt_time = time.time() - opt_start_time

        mip.write(f'{solstart_dir}/{lp_file_name}-start.sol')

        print(f'Obj: {mip.ObjVal:g}')
        
        obj_val = mip.ObjVal

    except gp.GurobiError as e:
        print(f'Error code {e.errno}: {e}')
        obj_val = None

    except AttributeError:
        print('Encountered an attribute error')
        obj_val = None

    if obj_val is not None:
        return obj_val, lp_file_name, opt_time
    else:
        print("Object value was not assigned due to an error.")

def optimize_multi_start(seq_files, seq_number, lpstart_dir, solstart_dir, s_start):
    
    chain_file = seq_files[seq_number]
    chain_name_with_ext = os.path.basename(chain_file)        
    chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
    lp_file_name = chain_name_without_ext
    seq_data = parse_seq_file(chain_file)
    this_RNA = seq_data['sequence']

    print(chain_name_without_ext)
    print(this_RNA)

    subseq_len = round(len(this_RNA)*0.63)

    mip = gp.Model(f'MIP-{seq_number}')

    # set the number of loops
    numH = 1#len(this_RNA)//5
    numI = 3#len(this_RNA)//4
    numB = 0#len(this_RNA)//3
    numM = 0#len(this_RNA)//3

    start = 1

    try:
        
        listP, listQ, listH, listI, listX = add_binary_vars(this_RNA, mip, start)
        mip.update()

        mip.setObjective(objectiveStartTerm(this_RNA, listQ, listH, listI), GRB.MINIMIZE)
        # mip.addConstr(objectiveStartTerm(this_RNA, listQ, listH, listI) >= MFE,"CMFE")
        # mip.addConstr(mip.getVarByName(f'B(14,15,28,30)') == 1)
        # mip.addConstr(mip.getVarByName(f'H(19,24)') == 1)
        # mip.addConstr(mip.getVarByName(f'I(6,12,37,41)') == 1)
        onePairConstraints(this_RNA, mip)
        noCrossConstraints(this_RNA, mip)
        stemConstraints(this_RNA, mip)
        # firstPairConstraints(this_RNA, mip)
        # lastPairConstraints(this_RNA, mip)
        # hairpinNTConstraints(this_RNA, mip)
        hairpinZeroConstraints(this_RNA, mip, maxH)
        hairpinIfThenConstraints(this_RNA, mip)
        # hairpinOnlyIfConstraints(this_RNA, mip)
        # numHairpinConstraints(this_RNA, numH, mip)
        # internalNTConstraints(this_RNA, mip)
        internalZeroConstraints(this_RNA, mip, maxI)
        internalIfThenConstraints(this_RNA, mip)
        internalOnlyIfConstraints(this_RNA, mip)
        # numInternalConstraints(this_RNA, numI, mip)
        # bulgeNTConstraints(this_RNA, mip)
        # bulgeZeroConstraints(this_RNA, mip, maxB)
        # bulgeIfThenConstraints(this_RNA, mip)
        # bulgeOnlyIfConstraints(this_RNA, mip)
        # numBulgeConstraints(this_RNA, numB, mip)
        # multiNTConstraints(this_RNA, mip)
        # multiZeroConstraints(this_RNA,mip,maxM)
        # multiIfThenConstraints(this_RNA, mip)
        # multiOnlyIfConstraints(this_RNA, mip)
        # numMultiConstraints(this_RNA, numM, mip)
        # numGreaterMultiConstraints(this_RNA, numM, mip)
        # consecutiveUnpairedConstraints(this_RNA, L, mip)
        loopNTConstraints(this_RNA, mip)        

        # add start constaints
        stemStartConstraints(this_RNA, mip, s_start, subseq_len+s_start)
        hairpinStartConstraints(this_RNA, mip, s_start, subseq_len+s_start)
        internalStartConstraints(this_RNA, mip, s_start, subseq_len+s_start)

        mip.update()   

        start_time = time.time()
        mip.write(f'{lpstart_dir}/{lp_file_name}-start-{s_start}.lp')
        print(f"--- {time.time() - start_time} seconds ---")

        # uncomment optimization tuning params if needed
        # mip.setParam("MIPGap", 0.2)
        # mip.setParam("PoolSolutions", 20)
        # mip.setParam("PoolSearchMode", 2)
        # mip.setParam("SolFiles", f"{chain_f}-decomposition-grb")
        mip.setParam("LogFile", f'{grb_log_dir}/log-{lp_file_name}-{s_start}')
        # mip.setParam("MIPFocus", 2)
        # mip.setParam("ConcurrentMIP ", 5)
        # mip.setParam("Presolve", 2)
        # mip.setParam("Cuts", 3)
        # mip.setParam("CliqueCuts", 2)

        # vars_with_Q = [var for var in mip.getVars() if 'Q' in var.VarName]
        # for i in range(len(vars_with_Q)):
        #     vars_with_Q[i].setAttr("BranchPriority",100)
        
        
        # mip.update()
        
        opt_start_time = time.time()        
        mip.optimize()
        opt_time = time.time() - opt_start_time

        mip.write(f'{solstart_dir}/{lp_file_name}-start-{s_start}.sol')

        print(f'Obj: {mip.ObjVal:g}')
        
        obj_val = mip.ObjVal

    except gp.GurobiError as e:
        print(f'Error code {e.errno}: {e}')
        obj_val = None

    except AttributeError:
        print('Encountered an attribute error')
        obj_val = None

    if obj_val is not None:
        return obj_val, lp_file_name, opt_time
    else:
        print("Object value was not assigned due to an error.")

#### TEST #######
# seq_number = 1
# chain_file = seq_files[seq_number]
# chain_name_with_ext = os.path.basename(chain_file)        
# chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
# lp_file_name = chain_name_without_ext
# seq_data = parse_seq_file(chain_file)
# RNA = seq_data['sequence']

# print(chain_name_without_ext)
# print(RNA)

# print(G_hairpin(RNA, 16,31))
# print(G_bulge(RNA, 8,9,38,43))


# mip = gp.Model(f'MIP-{seq_number}')

# add_binary_vars(RNA, mip)
# bulgeNTConstraints(RNA, mip)
# bulgeZeroConstraints(RNA, mip, maxB)
# bulgeIfThenConstraints(RNA, mip)
# bulgeOnlyIfConstraints(RNA, mip)
# numBulgeConstraints(RNA, noB, mip)



