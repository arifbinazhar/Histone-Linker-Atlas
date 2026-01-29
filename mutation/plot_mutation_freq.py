import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("data/processed/gbm/gbm_H1_mutation_summary.csv")

plt.figure(figsize=(8,5))
plt.bar(df["Gene"], df["Mutation_Count"])
plt.xticks(rotation=45)
plt.title("H1 Gene Mutation Frequency in Glioblastoma (TCGA)")
plt.ylabel("Mutation Count")
plt.tight_layout()
plt.savefig("data/processed/gbm/gbm_H1_mutation_barplot.png")
