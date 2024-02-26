from binary_variables import *
from constraints import *
from objective import *

OUT = open(f'{chain_f}-decomposition.lp', "w")

OUT.write("Maximize\n")
OUT.write("E \n")

OUT.write("\nsuch that \n")

objective = stemTerm(RNA) + fTerm(RNA) + lTerm(RNA) + hairpinTerm(RNA) + internalTerm(RNA) + bulgeTerm(RNA)

OUT.write(f'{objective} - E = 0 \n\n')

constraints = onePairConstraints(RNA) + noCrossConstraints(RNA) + stemConstraints(RNA) + firstPairConstraints(RNA) + lastPairConstraints(RNA) + hairpinOnlyIfConstraints(RNA) + hairpinIfThenConstraints(RNA) + hairpinNTConstraints(RNA)+ internalOnlyIfConstraints(RNA) + internalIfThenConstraints(RNA) + internalNTConstraints(RNA)  + bulgeOnlyIfConstraints(RNA) + bulgeIfThenConstraints(RNA) + bulgeNTConstraints(RNA) + numHairpinConstraints(RNA, numH) + numInternalConstraints(RNA, numI) + numBulgeConstraints(RNA, numB)

OUT.write(constraints)

OUT.write("\n")
OUT.write ("binary \n")
OUT.write("\n")

variables = listP + listQ + listF + listL + listH + listI + listX + listY + listB + listZ

OUT.write(variables)

OUT.write("\n")
OUT.write("general\n")
OUT.write("\n")
OUT.write("E\n")
OUT.write("\n")
OUT.write("end \n")
        
OUT.close()
