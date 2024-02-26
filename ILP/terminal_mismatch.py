import pandas as pd
import numpy as np

# Terminal mismatches are non-canonical pairs adjacent to helix ends.
c = 100

bps = np.array(['AU','CG','GC','GU','UA','UG'])

nucleotides = np.array(['A', 'C', 'G', 'U'])

mismatch = np.array([[[-0.8, -1.0, -0.8, -1.0],[-0.6, -0.7, -0.6, -0.7],[-0.8, -1.0, -0.8, -1.0],[-0.6, -0.8, -0.6, -0.8]],
                     [[-1.5, -1.5, -1.4, -1.5],[-1.0, -1.1, -1.0, -0.8],[-1.4, -1.5, -1.6, -1.5],[-1.0, -1.4, -1.0, -1.2]],
                     [[-1.1, -1.5, -1.3, -1.5],[-1.1, -0.7, -1.1, -0.5],[-1.6, -1.5, -1.4, -1.5],[-1.1, -1.0, -1.1, -0.7]],
                     [[-0.3, -1.0, -0.8, -1.0],[-0.6, -0.7, -0.6, -0.7],[-0.6, -1.0, -0.8, -1.0],[-0.6, -0.8, -0.6, -0.6]],
                     [[-1.0, -0.8, -1.1, -0.8],[-0.7, -0.6, -0.7, -0.5],[-1.1, -0.8, -1.2, -0.8],[-0.7, -0.6, -0.7, -0.5]],
                     [[-1.0, -0.8, -1.1, -0.8],[-0.7, -0.6, -0.7, -0.5],[-0.5, -0.8, -0.8, -0.8],[-0.7, -0.6, -0.7, -0.5]]])


mismatch_data = mismatch.reshape(bps.size*nucleotides.size, nucleotides.size) * c

# print(nucleotides)
# print(bps)
# print(mismatch_data)

midx = pd.MultiIndex.from_product([bps, nucleotides])

mismatch_df = pd.DataFrame(mismatch_data, index = midx, columns=nucleotides) 

# print(mismatch_df)

# print(mismatch_df.loc['AU', 'C']['A'])







