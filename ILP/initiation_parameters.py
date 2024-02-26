import pandas as pd
import numpy as np

# DESTABILIZING ENERGIES BY SIZE OF LOOP (INTERPOLATE WHERE NEEDED)

loops = ['internal', 'bulge', 'hairpin']
n = np.arange(1,31)
c = 100

# print(n)

energies = np.array([[ 0, 3.8, 0],
                [ 0, 2.8, 0],
                [ 0, 3.2, 5.4],
                [1.1, 3.6, 5.6],
                [2.0, 4.0, 5.7],
                [2.0, 4.4, 5.4],
                [2.1, 4.6, 6.0],
                [2.3, 4.7, 5.5],
                [2.4, 4.8, 6.4],
                [2.5, 4.9, 6.5],
                [2.6, 5.0, 6.6],
                [2.7, 5.1, 6.7],
                [2.8, 5.2, 6.8],
                [2.9, 5.3, 6.9],
                [2.9, 5.4, 6.9],
                [3.0, 5.4, 7.0],
                [3.1, 5.5, 7.1],
                [3.1, 5.5, 7.1],
                [3.2, 5.6, 7.2],
                [3.3, 5.7, 7.2],
                [3.3, 5.7, 7.3],
                [3.4, 5.8, 7.3],
                [3.4, 5.8, 7.4],
                [3.5, 5.8, 7.4],
                [3.5, 5.9, 7.5],
                [3.5, 5.9, 7.5],
                [3.6, 6.0, 7.5],
                [3.6, 6.0, 7.6],
                [3.7, 6.0, 7.6],
                [3.7, 6.1, 7.7]])

dados = energies * c
linhas = np.arange(1,31)
colunas = loops

initiation_df = pd.DataFrame(
    data=dados,
    index=linhas,
    columns=colunas
)

# print('\nDataFrame:', initiation_df)

# print(initiation_df.loc[ 3, 'hairpin'])