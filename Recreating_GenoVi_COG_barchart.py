import pandas as pd
import matplotlib.pyplot as plt
plt.style.use("seaborn-v0_8-darkgrid")

# Data entry
Counts = [
    ("D", 111), ("M", 307), ("N", 57), ("O", 227), ("T", 382), ("U", 95), ("V", 184), ("W", 51), ("Y", 0), ("Z", 3), ("A", 0), ("B", 0), ("J", 243), ("K", 280), ("L", 149), ("X", 153), ("C", 328), ("E", 209), ("F", 79), ("G", 178), ("H", 248), ("I", 145), ("P", 199), ("Q", 91),("R", 471), ("S", 253), ("Unclassified", 20)
]

df = pd.DataFrame(Counts, columns=["COG", "Count"])

# Replacing the characters with the COG classificiation
cog_names = {
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "X": "Mobilome: prophages, transposons",
    "C": "Energy production and conversion",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "Unclassified": "Unclassified"
}


df["Label"] = df["COG"].map(cog_names)

# Strong color palette missing
color_palette = [
 "#78eded"
]

df["Color"] = [color_palette[i % len(color_palette)] for i in range(len(df))]

# Building of the plot
fig, ax = plt.subplots(figsize=(14, 6), constrained_layout=False)

# X-axis Labelling
bars = ax.bar(
    df["Label"],
    df["Count"],
    color=df["Color"],
    width=0.8
)
ax.set_title("COG Category", fontsize=14)
ax.set_xticks(range(len(df)))
ax.set_xticklabels(
    df["Label"],
    rotation=65,
    ha='right',
    va='top',
    rotation_mode='anchor'
)
ax.tick_params(axis="x", labelsize=8, pad=10)
ax.margins(x=0.02)
plt.subplots_adjust(bottom=0.36)

# Values on bars
for bar in bars:
    height = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        height + 5,
        f"{int(height)}",
        ha="center",
        va="bottom",
        fontsize=8
    )

plt.tight_layout()
plt.show()





