from binary_variables_grb import *
from constraints_grb import *
from objective_grb import *


try:
    mip.setObjective(objectiveTerm(RNA), GRB.MINIMIZE)

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

    mip.update()

    mip.write(f'{chain_f}-decomposition-grb.lp')

    mip.optimize()

    mip.write(f'{chain_f}-decomposition-grb.sol')

    print(f'Obj: {mip.ObjVal:g}')

except gp.GurobiError as e:
    print(f'Error code {e.errno}: {e}')

except AttributeError:
    print('Encountered an attribute error')


# OUT = open(f'{chain_f}-decomposition.lp', "w")

# OUT.write("Maximize\n")
# OUT.write("E \n")

# OUT.write("\nsuch that \n")

# objective = stemTerm(RNA) + hairpinTerm(RNA) + internalTerm(RNA) + bulgeTerm(RNA) #+ fTerm(RNA) + lTerm(RNA)

# OUT.write(f'{objective} - E = 0 \n\n')

# constraints = onePairConstraints(RNA) + noCrossConstraints(RNA) + stemConstraints(RNA) + firstPairConstraints(RNA) + lastPairConstraints(RNA) + hairpinOnlyIfConstraints(RNA) + hairpinIfThenConstraints(RNA) + hairpinNTConstraints(RNA)+ internalOnlyIfConstraints(RNA) + internalIfThenConstraints(RNA) + internalNTConstraints(RNA)  + bulgeOnlyIfConstraints(RNA) + bulgeIfThenConstraints(RNA) + bulgeNTConstraints(RNA) + numHairpinConstraints(RNA, numH) + numInternalConstraints(RNA, numI) + numBulgeConstraints(RNA, numB)

# OUT.write(constraints)

# OUT.write("\n")
# OUT.write ("binary \n")
# OUT.write("\n")

# variables = listP + listQ + listF + listL + listH + listI + listX + listY + listB + listZ

# OUT.write(variables)

# OUT.write("\n")
# OUT.write("general\n")
# OUT.write("\n")
# OUT.write("E\n")
# OUT.write("\n")
# OUT.write("end \n")
        
# OUT.close()
