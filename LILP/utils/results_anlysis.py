import pandas as pd
from constants_paths import *

model_name = 'lilp-nomulti'
results_file = f'{results_folder_name}/LILP_{len_start}_{len_end}_{model_name}.txt'
df = pd.read_csv(f'{results_file}', sep=r"\s+")

print(df.head())

import seaborn as sns
import matplotlib.pyplot as plt

# Example: violin plot for F1 and Fbeta
data=df[["INFILP", "INFRNAstr", "INFRNAFold", "INFUNAFold"]]
means = data.mean()

sns.violinplot(data=data, inner=None)

# for i, m in enumerate(means):
#     plt.scatter(i, m, color="red", s=40, zorder=3)

for i, m in enumerate(means):
    plt.text(i, m, f"{m:.3f}", color="black", ha="center", va="bottom", fontweight="bold")


plt.title("Distribution of F1 and Fbeta Scores")
plt.ylabel("Score Value")
plt.ylim(0, 1.0)
plt.show()


# df_melted = df.melt(id_vars=["RNAseqname"], value_vars=["F1ILP", "F1RNAstr", "F1UNAFold"],
#                     var_name="Metric", value_name="Value")

# sns.violinplot(x="Metric", y="Value", data=df_melted)
# plt.title("RNA Structure Prediction Scores")
# plt.ylim(0, 1.0)
# plt.show()
