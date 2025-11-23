import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from constants_paths import *

model_name = 'lilp-multi'
# results_nomulti_file = f'{results_folder_name}/LILP_50_80_nomulti.txt'
# results_multi_file = f'{results_folder_name}/LILP_50_80_multi.txt'
results_multi_file = f'{results_folder_name}/LILP_70_80_{model_name}.txt'
# f'{results_folder_name}/LILP_{len_start}_{len_end}_{model_name}.txt'
# df_nomulti = pd.read_csv(f'{results_nomulti_file}', sep=r"\s+")
df_multi = pd.read_csv(f'{results_multi_file}', sep=r"\s+")

#print(df_nomulti)

# cols_nomulti = ["INFILP"]
cols_multi = ["INFILP", "INFRNAstr", "INFRNAFold", "INFUNAFold"]

# data_nomulti = df_nomulti[cols_nomulti]
data_multi = df_multi[cols_multi]

# data_nomulti = data_nomulti.rename(columns={'INFILP': '2-degree LILP'})
data_multi = data_multi.rename(columns={'INFILP': '3-degree LILP', 'INFRNAstr': 'RNAstucture', 'INFRNAFold' : 'RNAFold', 'INFUNAFold': 'UNAFold'})

# combined = pd.concat([data_nomulti, data_multi], axis=1)
combined = data_multi

means = combined.mean()

plt.figure(figsize=(8,8))
sns.violinplot(data=combined, inner='point')
for i, m in enumerate(means):
    plt.text(i, m, f"{m:.3f}", color="black", ha="center", va="bottom", fontweight="bold")
plt.title("Distribution of INF metric")
plt.ylabel("INF Value")
plt.xticks(rotation=30)
plt.ylim(0, 1.0)
plt.tight_layout()
plt.show()

# df_melted = df.melt(id_vars=["RNAseqname"], value_vars=["F1ILP", "F1RNAstr", "F1UNAFold"],
#                     var_name="Metric", value_name="Value")

# sns.violinplot(x="Metric", y="Value", data=df_melted)
# plt.title("RNA Structure Prediction Scores")
# plt.ylim(0, 1.0)
# plt.show()
