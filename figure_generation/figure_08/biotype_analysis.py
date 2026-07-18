
"""
Preprocessing steps starting with unaligned bam files after basecalling

Alignment to reference:
(For rRNA and polyA)
samtools fastq -T "*" in.ubam | minimap2 -y --MD -t 5 -ax splice-ont -uf -k14 /global/cfs/cdirs/m5243/references/hg38.fa -  | samtools sort -@ 2 | samtools view --threads 2 -hb -F 260 > out.aligned.sorted.bam

(For tRNA)
samtools fastq -T "*" in.ubam | minimap2 -y -ax splice -t 5 -k 7 -w 3 -A 2 -B 1 -O1,32 -E2,0 -n 1 -m 13 -s 30 --secondary=no --MD /global/cfs/cdirs/m5243/references/hg38.rtRNA.with_oligos.fa - | samtools sort -@ 5 | samtools view --threads 5 -hb > out.rtRNA.all.sorted.bam
samtools -hbS -F 4 out.rtRNA.all.sorted.bam > out.rtRNA.aligned.sorted.bam (Collect rRNA and tRNA alinging reads)
samtools -hbS -f 4 out.rtRNA.all.sorted.bam > out.rtRNA.unaligned.sorted.bam (Collect unaligned reads)
samtools fastq -T "*" out.rtRNA.unaligned.sorted.bam | minimap2 -y --MD -t 5 -ax splice-ont -uf -k14 /global/cfs/cdirs/m5243/references/hg38.fa -  | samtools sort -@ 2 | samtools view --threads 2 -hb -F 260 > out.genome.aligned.bam
samtools merge -@ 5 out.aligned.sorted.bam out.genome.aligned.bam > out.aligned.sorted.bam

Quantification with feature counts:
featureCounts -a combined_v49_rRNAs.gtf -L -T 50 -o "out.aligned.sorted.tsv" out.aligned.sorted.bam
"""

#%%

"""
polyA RNA
"""
import pandas as pd
import polars as pl
import dmode
import matplotlib.pyplot as plt


gene_type_colors = {
    # Protein-coding — deep forest greens
    "protein_coding": "#2E5E3A",
    "TEC":            "#6B8E5A",

    # tRNA family — warm earth/bark tones
    "tRNA":           "#A0522D",
    "Mt_tRNA":        "#8B4513",

    # rRNA family — clay / terracotta
    "rRNA":           "#C1693C",
    "Mt_rRNA":        "#9C4A2A",
    "rRNA_pseudogene":"#D4926A",
    "18S":            "#E08D5B",
    "28S":            "#C76B3A",
    "5.8S":           "#E8A87C",
    "5S":             "#B5651D",

    # Small/regulatory RNAs — wildflower & moss tones
    "miRNA":          "#7B68A6",   # lavender
    "snoRNA":         "#9370B5",   # mauve-violet
    "snRNA":          "#5D8AA8",   # cornflower
    "scaRNA":         "#6CA0B5",   # dusty teal
    "sRNA":           "#88B04B",   # moss green
    "misc_RNA":       "#A3A847",   # olive
    "vault_RNA":      "#4F7942",   # fern
    "ribozyme":       "#3D6B4A",   # pine
    "lncRNA":         "#4A7C8C",   # slate teal

    # IG genes — autumn reds/oranges
    "IG_C_gene":      "#B22222",
    "IG_V_gene":      "#C1432E",
    "IG_D_gene":      "#D35438",
    "IG_J_gene":      "#A83232",

    # IG pseudogenes — faded autumn
    "IG_C_pseudogene":"#D98C7A",
    "IG_V_pseudogene":"#E0A08F",
    "IG_J_pseudogene":"#CD9384",
    "IG_pseudogene":  "#C99182",

    # TR genes — golden / amber
    "TR_C_gene":      "#DAA520",
    "TR_V_gene":      "#E0B23A",
    "TR_D_gene":      "#C99700",
    "TR_J_gene":      "#EAC54F",

    # TR pseudogenes — muted gold
    "TR_V_pseudogene":"#D4BC72",
    "TR_J_pseudogene":"#CBB76A",

    # General pseudogenes — stone / lichen greys
    "processed_pseudogene":              "#9A9B7A",
    "unprocessed_pseudogene":            "#8A8B6C",
    "transcribed_processed_pseudogene":  "#A8A98A",
    "transcribed_unprocessed_pseudogene":"#7C7D60",
    "transcribed_unitary_pseudogene":    "#B0B190",
    "unitary_pseudogene":                "#90916F",
    "translated_processed_pseudogene":   "#A0A17E",

    # Other 
    "artifact":       "#6E6E6E",   # neutral grey
}


#%%

base = "/global/cfs/cdirs/m5243/final_bedRmods/ont/read_count_tables/"

HRP_A_001_001_df = pd.read_csv(base + "HRP_A_001_native_polyA_RNA_001.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")

HRP_A_002_001_df = pd.read_csv(base + "HRP_A_002_native_polyA_RNA_001.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_003_001_df = pd.read_csv(base + "HRP_A_003_native_polyA_RNA_001.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_005_001_df = pd.read_csv(base + "HRP_A_005_native_polyA_RNA_001.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_006_004_df = pd.read_csv(base + "HRP_A_006_native_polyA_RNA_004.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_007_001_df = pd.read_csv(base + "HRP_A_007_native_polyA_RNA_001.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_009_001_df = pd.read_csv(base + "HRP_A_009_native_polyA_RNA_001.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_009_002_df = pd.read_csv(base + "HRP_A_009_native_polyA_RNA_002.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_010_001_df = pd.read_csv(base + "HRP_A_010_native_polyA_RNA_001.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_010_002_df = pd.read_csv(base + "HRP_A_010_native_polyA_RNA_002.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_011_001_df = pd.read_csv(base + "HRP_A_011_native_polyA_RNA_001.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_011_002_df = pd.read_csv(base + "HRP_A_011_native_polyA_RNA_002.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_011_003_df = pd.read_csv(base + "HRP_A_011_native_polyA_RNA_003.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_013_001_df = pd.read_csv(base + "HRP_A_013_native_polyA_RNA_001.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_014_001_df = pd.read_csv(base + "HRP_A_014_native_polyA_RNA_001.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
combined_df = pd.read_csv(base + "polyA_native_merged.merged.tsv", header=0, sep="\t" ,comment="#")

dfs = {
    "HRP_A_001_1": HRP_A_001_001_df,
    "HRP_A_002_1": HRP_A_002_001_df,
    "HRP_A_003_1": HRP_A_003_001_df,
    "HRP_A_005_1": HRP_A_005_001_df,
    "HRP_A_006_4": HRP_A_006_004_df,
    "HRP_A_007_1": HRP_A_007_001_df,
    "HRP_A_009_1": HRP_A_009_001_df,
    "HRP_A_009_2": HRP_A_009_002_df,
    "HRP_A_010_1": HRP_A_010_001_df,
    "HRP_A_010_2": HRP_A_010_002_df,
    "HRP_A_011_1": HRP_A_011_001_df,
    "HRP_A_011_2": HRP_A_011_002_df,
    "HRP_A_011_3": HRP_A_011_003_df,
    "HRP_A_013_1": HRP_A_013_001_df,
    "HRP_A_014_1": HRP_A_014_001_df,
    "combined": combined_df
}

for name, df in dfs.items():
    df.rename(columns={df.columns[-1]: name}, inplace=True)
merged = pd.concat(
    [df.set_index("Geneid")[name] for name, df in dfs.items()],
    axis=1,
)

gtf_df, gene_body_df = dmode.metagene_plot.prepare_gene_body_coverage("/global/cfs/cdirs/m5243/references/gencode.v49.basic.annotation.gtf")

# %%
# 1. reads per kilobase
#rpk = merged.div(length / 1000, axis=0)
# 2. per-sample scaling factor (sum of RPK, in millions)
scaling = merged.sum(axis=0) / 1e6
# 3. divide each column by its scaling factor
cpm = merged.div(scaling, axis=1)

#%%
gtf_df_selection = gtf_df.filter(pl.col("feature") == "gene").select("gene_type", "gene_name", "gene_id").to_pandas()
merged_counts_and_gtf_df =  pd.merge(cpm,gtf_df_selection,left_on=cpm.index, right_on=gtf_df_selection["gene_id"], how="inner")

#%%
pd.set_option('display.float_format', '{:.6f}'.format)
source_table = merged_counts_and_gtf_df.groupby("gene_type").sum().drop(columns=["gene_id", "gene_name", "key_0"]).sort_values(by="combined", ascending=False)
source_table = (source_table / 1000000) * 100
source_table.to_csv(base + "gene_type_composition_polyA_combined_only.tsv", sep="\t",index=True, header=True)

#%%

plt.rcParams.update({
    "font.size": 20,          # base
    "axes.titlesize": 24,
    "axes.labelsize": 28,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 20,
    "legend.title_fontsize": 22,
})

#%%
sample_cols = ["combined"
]

comp = merged_counts_and_gtf_df.groupby("gene_type")[sample_cols].sum()   # absolute cpm
rel = comp.div(comp.sum(axis=0), axis=1)                # relative widths

# keep gene_types exceeding 2% relative share in at least one sample
keep = rel.index[(rel > 0.01).any(axis=1)]

# order by total signal (so "other" lands last and bars stack consistently)
keep = rel.loc[keep].sum(axis=1).sort_values(ascending=False).index

def collapse(df):
    out = df.loc[keep].copy()
    out.loc["other"] = df.drop(keep).sum(axis=0)
    return out

rel_plot = collapse(rel)      # controls segment width (0–1)
abs_plot = collapse(comp)     # absolute cpm for labels


fig,ax = plt.subplots()
colors = [gene_type_colors.get(b, "#CCCCCC") for b in rel_plot.index]

rel_plot.T.plot(
    kind="barh", stacked=True, figsize=(16,4),
    color=colors, width=0.8, ax=ax
)


ax.set_xlabel("Relative composition (CPM / 1M)")
ax.set_ylabel("polyA RNA")
ax.set_yticklabels("")

ax.set_xlim(0, 1)
ax.legend(title="Biotype", bbox_to_anchor=(1.01, 1), loc="upper left")


plt.tight_layout()

plt.savefig(base + "Figure1_gene_type_composition_combined_only_CPMS.pdf", format="pdf", bbox_inches="tight")

plt.show()

# %%
sample_cols = [
    "HRP_A_001_1","HRP_A_002_1","HRP_A_003_1","HRP_A_005_1",
    "HRP_A_006_4","HRP_A_007_1","HRP_A_009_1","HRP_A_009_2",
    "HRP_A_010_1","HRP_A_010_2","HRP_A_011_1","HRP_A_011_2",
    "HRP_A_011_3","HRP_A_013_1","HRP_A_014_1"
][::-1]

comp = merged_counts_and_gtf_df.groupby("gene_type")[sample_cols].sum()   # absolute cpm
rel = comp.div(comp.sum(axis=0), axis=1)                # relative widths

# keep gene_types exceeding 2% relative share in at least one sample
keep = rel.index[(rel > 0.01).any(axis=1)]

# order by total signal (so "other" lands last and bars stack consistently)
keep = rel.loc[keep].sum(axis=1).sort_values(ascending=False).index

def collapse(df):
    out = df.loc[keep].copy()
    out.loc["other"] = df.drop(keep).sum(axis=0)
    return out


rel_plot = collapse(rel)      # controls segment width (0–1)
abs_plot = collapse(comp)     # absolute cpm for labels

fig,ax = plt.subplots()
colors = [gene_type_colors.get(b, "#CCCCCC") for b in rel_plot.index]

rel_plot.T.plot(
    kind="barh", stacked=True, figsize=(16,8),
    color=colors, width=0.8, ax=ax
)

ax.set_xlabel("Relative composition (CPM / 1M)")
ax.set_ylabel("Sample")
ax.set_xlim(0, 1)
ax.legend(title="Biotype", bbox_to_anchor=(1.01, 1), loc="upper left")

plt.savefig("Figure1_gene_type_composition.png", dpi=300, bbox_inches="tight")
plt.savefig("Figure1_gene_type_composition.svg", format="svg", bbox_inches="tight")
plt.savefig("Figure1_gene_type_composition_CPMS.pdf", format="pdf", bbox_inches="tight")
plt.tight_layout()
plt.show()



"""
rRNA
"""
#%%
base = "/global/cfs/cdirs/m5243/analysis/unaligned_rRNA_bams/"
HRP_A_018_1_df = pd.read_csv(base + "native_rRNA_Novoa_Eppendorf_Rep1.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_018_2_df = pd.read_csv(base + "native_rRNA_Novoa_Eppendorf_Rep2.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_018_3_df = pd.read_csv(base + "native_rRNA_Novoa_Eppendorf_Rep3.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_017_1_df = pd.read_csv(base + "native_rRNA_Novoa_rRNA_native_rep1.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_017_2_df = pd.read_csv(base + "native_rRNA_Novoa_rRNA_native_rep2.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_017_3_df = pd.read_csv(base + "native_rRNA_Novoa_rRNA_native_rep3.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")
HRP_A_006_2_df = pd.read_csv(base + "native_rRNA_SUP_Gerber_native_rRNA_251118_r1_v1.merged.aligned.sorted.tsv",sep = "\t", header=0, comment="#")

combined_df = pd.read_csv(base + "merged_rRNA_samples.tsv",sep = "\t", header=0, comment="#")
# %%

# %%

# %%
dfs = {
    "HRP_A_006_2": HRP_A_006_2_df,
    "HRP_A_017_1": HRP_A_017_1_df,
    "HRP_A_017_2": HRP_A_017_2_df,
    "HRP_A_017_3": HRP_A_017_3_df,
    "HRP_A_018_1": HRP_A_018_1_df,
    "HRP_A_018_2": HRP_A_018_2_df,
    "HRP_A_018_3": HRP_A_018_3_df,
    "combined": combined_df
}

for name, df in dfs.items():
    df.rename(columns={df.columns[-1]: name}, inplace=True)
merged = pd.concat(
    [df.set_index("Geneid")[name] for name, df in dfs.items()],
    axis=1,
)




# %%
gtf_df, gene_body_df = dmode.metagene_plot.prepare_gene_body_coverage("/home/stefan/Synology/Data_jamboree/final_bedRmods/ont/unaligned_rRNA_bams/combined_v49_rRNAs.gtf")

# %%
print(merged)

# %%
# 1. reads per kilobase
#rpk = merged.div(length / 1000, axis=0)
# 2. per-sample scaling factor (sum of RPK, in millions)
scaling = merged.sum(axis=0) / 1e6
# 3. divide each column by its scaling factor
cpm = merged.div(scaling, axis=1)


#%%
gtf_df_selection = gtf_df.filter(pl.col("feature") == "gene").select("gene_type", "gene_name", "gene_id").to_pandas()

# %%
merged_counts_and_gtf_df =  pd.merge(cpm,gtf_df_selection,left_on=cpm.index, right_on=gtf_df_selection["gene_id"], how="inner")

merged_counts_and_gtf_df["gene_type"] = merged_counts_and_gtf_df["gene_type"].replace({
    "rRNA_5S": "rRNA",
    "rRNA_5.8S": "rRNA",
    "rRNA_18S": "rRNA",
    "rRNA_28S": "rRNA",
})

print(merged_counts_and_gtf_df)
#%%
print(merged_counts_and_gtf_df["gene_type"])

#%%
pd.set_option('display.float_format', '{:.6f}'.format)
source_table = merged_counts_and_gtf_df.groupby("gene_type").sum().drop(columns=["gene_id", "gene_name", "key_0"]).sort_values(by="combined", ascending=False)
source_table = (source_table / 1000000) * 100
source_table.to_csv(base + "gene_type_composition_rRNA_combined_only.tsv", sep="\t", index=True, header=True)

#%%
plt.rcParams.update({
    "font.size": 20,          # base
    "axes.titlesize": 24,
    "axes.labelsize": 28,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 18,
    "legend.title_fontsize": 22,
})
#%%
sample_cols = ["combined"]


comp = merged_counts_and_gtf_df.groupby("gene_type")[sample_cols].sum()   # absolute cpm
rel = comp.div(comp.sum(axis=0), axis=1)                # relative widths

# keep gene_types exceeding 2% relative share in at least one sample
keep = rel.index[(rel > 0.01).any(axis=1)]

# order by total signal (so "other" lands last and bars stack consistently)
keep = rel.loc[keep].sum(axis=1).sort_values(ascending=False).index

def collapse(df):
    out = df.loc[keep].copy()
    out.loc["other"] = df.drop(keep).sum(axis=0)
    return out

rel_plot = collapse(rel)      # controls segment width (0–1)
abs_plot = collapse(comp)     # absolute cpm for labels


fig,ax = plt.subplots()
colors = [gene_type_colors.get(b, "#CCCCCC") for b in rel_plot.index]

rel_plot.T.plot(
    kind="barh", stacked=True, figsize=(16,4),
    color=colors, width=0.8, ax=ax
)


ax.set_xlabel("Relative composition (CPM / 1M)")
ax.set_ylabel("rRNA")
ax.set_yticklabels("")

ax.set_xlim(0, 1)
ax.legend(title="Biotype", bbox_to_anchor=(1.01, 1), loc="upper left")


plt.tight_layout()
plt.show()
fig.savefig(base + "Figure1_rRNA_biotype_CPMS_single_dataset.pdf", format="pdf")


# %%
sample_cols = [
    "HRP_A_006_2",
    "HRP_A_017_1",
    "HRP_A_017_2",
    "HRP_A_017_3",
    "HRP_A_018_1",
    "HRP_A_018_2",
    "HRP_A_018_3"
][::-1]

comp = merged_counts_and_gtf_df.groupby("gene_type")[sample_cols].sum()   # absolute cpm
rel = comp.div(comp.sum(axis=0), axis=1)                # relative widths

# keep gene_types exceeding 2% relative share in at least one sample
keep = rel.index[(rel > 0.01).any(axis=1)]

# order by total signal (so "other" lands last and bars stack consistently)
keep = rel.loc[keep].sum(axis=1).sort_values(ascending=False).index

def collapse(df):
    out = df.loc[keep].copy()
    out.loc["other"] = df.drop(keep).sum(axis=0)
    return out


rel_plot = collapse(rel)      # controls segment width (0–1)
abs_plot = collapse(comp)     # absolute cpm for labels

fig,ax = plt.subplots()
colors = [gene_type_colors.get(b, "#CCCCCC") for b in rel_plot.index]

rel_plot.T.plot(
    kind="barh", stacked=True, figsize=(16,8),
    color=colors, width=0.8, ax=ax
)

ax.set_xlabel("Relative composition (CPM / 1M)")
ax.set_ylabel("Sample")
ax.set_xlim(0, 1)
ax.legend(title="Biotype", bbox_to_anchor=(1.01, 1), loc="upper left")
fig.savefig(base + "Figure1_rRNA_biotype_CPMS.pdf", format="pdf")

plt.tight_layout()
plt.show()



"""
tRNA
"""
#%%
base = "/global/cfs/cdirs/m5243/analysis/unaligned_rRNA_bams/realigned_bams/"
HRP_A_020_1_1_df = pd.read_csv(base + "HRP_A_020_1_native_tRNA_001_final_merge.tsv",sep = "\t", header=0, comment="#")
HRP_A_020_1_2_df = pd.read_csv(base + "HRP_A_020_1_native_tRNA_002_final_merge.tsv",sep = "\t", header=0, comment="#")
HRP_A_020_1_3_df = pd.read_csv(base + "HRP_A_020_1_native_tRNA_003_final_merge.tsv",sep = "\t", header=0, comment="#")
HRP_A_020_2_1_df = pd.read_csv(base + "HRP_A_020_2_native_tRNA_001_final_merge.tsv",sep = "\t", header=0, comment="#")
HRP_A_020_2_2_df = pd.read_csv(base + "HRP_A_020_2_native_tRNA_002_final_merge.tsv",sep = "\t", header=0, comment="#")
HRP_A_020_2_3_df = pd.read_csv(base + "HRP_A_020_2_native_tRNA_003_final_merge.tsv",sep = "\t", header=0, comment="#")
HRP_A_020_3_1_df = pd.read_csv(base + "HRP_A_020_3_native_tRNA_001_final_merge.tsv",sep = "\t", header=0, comment="#")
HRP_A_020_3_2_df = pd.read_csv(base + "HRP_A_020_3_native_tRNA_002_final_merge.tsv",sep = "\t", header=0, comment="#")
HRP_A_020_3_3_df = pd.read_csv(base + "HRP_A_020_3_native_tRNA_003_final_merge.tsv",sep = "\t", header=0, comment="#")
HRP_A_021_1_1_df = pd.read_csv(base + "HRP_A_021_1_native_tRNA_001_final_merge.tsv",sep = "\t", header=0, comment="#")
HRP_A_021_1_2_df = pd.read_csv(base + "HRP_A_021_1_native_tRNA_002_final_merge.tsv",sep = "\t", header=0, comment="#")
HRP_A_021_1_3_df = pd.read_csv(base + "HRP_A_021_1_native_tRNA_003_final_merge.tsv",sep = "\t", header=0, comment="#")
HRP_A_004_1_df = pd.read_csv(base + "native_tRNA_SUP_ADAMCZYK_M_NATIVE_tRNA_003_final_merge.tsv",sep = "\t", header=0, comment="#")
HRP_A_006_1_df = pd.read_csv(base + "native_tRNA_SUP_GERBER_ONT_native_tRNA_final_merge.tsv",sep = "\t", header=0, comment="#")
HRP_A_012_1_df = pd.read_csv(base + "native_tRNA_SUP_Soares_ONT_native_tRNA_final_merge.tsv",sep = "\t", header=0, comment="#")


combined_df = pd.read_csv(base + "/realigned_bams/combined_dataset_final_merge.tsv",sep = "\t", header=0, comment="#")

# %%
dfs = {
"HRP_A_004_1" : HRP_A_004_1_df,
"HRP_A_006_1" : HRP_A_006_1_df,
"HRP_A_012_1" : HRP_A_012_1_df,
"HRP_A_020_1_1" : HRP_A_020_1_1_df,
"HRP_A_020_1_2" : HRP_A_020_1_2_df,
"HRP_A_020_1_3" : HRP_A_020_1_3_df,
"HRP_A_020_2_1" : HRP_A_020_2_1_df,
"HRP_A_020_2_2" : HRP_A_020_2_2_df,
"HRP_A_020_2_3" : HRP_A_020_2_3_df,
"HRP_A_020_3_1" : HRP_A_020_3_1_df,
"HRP_A_020_3_2" : HRP_A_020_3_2_df,
"HRP_A_020_3_3" : HRP_A_020_3_3_df,
"HRP_A_021_1_1" : HRP_A_021_1_1_df,
"HRP_A_021_1_2" : HRP_A_021_1_2_df,
"HRP_A_021_1_3" : HRP_A_021_1_3_df,
"combined" :combined_df
}

for name, df in dfs.items():
    df.rename(columns={df.columns[-1]: name}, inplace=True)
merged = pd.concat(
    [df.set_index("Geneid")[name] for name, df in dfs.items()],
    axis=1,
)

# %%
gtf_df, gene_body_df = dmode.metagene_plot.prepare_gene_body_coverage("/mnt/data01/stpastore/final_bedRmods/unaligned_tRNA_bams/combined_v49_rRNAs_tRNAs.gtf")



# %%
# 1. reads per kilobase
#rpk = merged.div(length / 1000, axis=0)
# 2. per-sample scaling factor (sum of RPK, in millions)
scaling = merged.sum(axis=0) / 1e6
# 3. divide each column by its scaling factor
cpm = merged.div(scaling, axis=1)


#%%
print(cpm)

#%%
print(gtf_df.columns)
#%%
import polars as pl
gtf_df_selection = gtf_df.filter(pl.col("feature") == "gene").select("gene_type", "gene_name", "gene_id").to_pandas()

# %%
merged_counts_and_gtf_df =  pd.merge(cpm,gtf_df_selection,left_on=cpm.index, right_on=gtf_df_selection["gene_id"], how="inner")

merged_counts_and_gtf_df["gene_type"] = merged_counts_and_gtf_df["gene_type"].replace({
    # rRNA
    "rRNA_5S": "5S",
    "rRNA_5.8S": "5.8S",
    "rRNA_18S": "18S",
    "rRNA_28S": "28S",
    "mt-rRNA_12s": "mt-12S",
    "mt-rRNA_16s": "mt-16S",
    # cytoplasmic tRNA
    "tRNAAla_CGC": "Ala_CGC",
    "tRNAAla_TGC": "Ala_TGC",
    "tRNAAla_AGC": "Ala_AGC",
    "tRNAArg_ACG": "Arg_ACG",
    "tRNAArg_TCG": "Arg_TCG",
    "tRNAArg_CCG2": "Arg_CCG2",
    "tRNAArg_TCT2": "Arg_TCT2",
    "tRNAArg_CCT": "Arg_CCT",
    "tRNAArg_CCG1": "Arg_CCG1",
    "tRNAArg_TCT1": "Arg_TCT1",
    "tRNAAsn_GTT": "Asn_GTT",
    "tRNAAsp_GTC": "Asp_GTC",
    "tRNACys_GCA": "Cys_GCA",
    "tRNAGln_CTG_TTG": "Gln_CTG_TTG",
    "tRNAGlu_CTC": "Glu_CTC",
    "tRNAGlu_TTC": "Glu_TTC",
    "tRNAGly_GCC": "Gly_GCC",
    "tRNAGly_TCC": "Gly_TCC",
    "tRNAGly_CCC": "Gly_CCC",
    "tRNAHis_GTG": "His_GTG",
    "tRNAIle_AAT": "Ile_AAT",
    "tRNAIle_TAT": "Ile_TAT",
    "tRNAIle_GAT": "Ile_GAT",
    "tRNALeu_AAG": "Leu_AAG",
    "tRNALeu_CAG": "Leu_CAG",
    "tRNALeu_CAA": "Leu_CAA",
    "tRNALeu_TAA": "Leu_TAA",
    "tRNALys_CTT": "Lys_CTT",
    "tRNALys_TTT": "Lys_TTT",
    "tRNAMet_CAT": "Met_CAT",
    "tRNAPhe_GAA": "Phe_GAA",
    "tRNAPro_AGG_CGG_TGG": "Pro_AGG_CGG_TGG",
    "tRNASeC_TCA": "SeC_TCA",
    "tRNASer_AGA": "Ser_AGA",
    "tRNASer_CGA": "Ser_CGA",
    "tRNASer_GCT": "Ser_GCT",
    "tRNAThr_AGT": "Thr_AGT",
    "tRNAThr_CGT1": "Thr_CGT1",
    "tRNAThr_TGT": "Thr_TGT",
    "tRNAThr_CGT2": "Thr_CGT2",
    "tRNATrp_CCA": "Trp_CCA",
    "tRNATyr_GTA2": "Tyr_GTA2",
    "tRNATyr_ATA": "Tyr_ATA",
    "tRNATyr_GTA1": "Tyr_GTA1",
    "tRNAVal_AAC_CAC": "Val_AAC_CAC",
    "tRNAVal_TAC": "Val_TAC",
    # mitochondrial tRNA
    "mttRNAAla_TGC": "mt-Ala_TGC",
    "mttRNAArg_TCG": "mt-Arg_TCG",
    "mttRNAAsn_GTT": "mt-Asn_GTT",
    "mttRNAAsp_GTC": "mt-Asp_GTC",
    "mttRNACys_GCA": "mt-Cys_GCA",
    "mttRNAGln_TTG": "mt-Gln_TTG",
    "mttRNAGlu_TTC": "mt-Glu_TTC",
    "mttRNAGly_TCC": "mt-Gly_TCC",
    "mttRNAHis_GTG": "mt-His_GTG",
    "mttRNAIle_GAT": "mt-Ile_GAT",
    "mttRNALeu_TAA": "mt-Leu_TAA",
    "mttRNALeu_TAG": "mt-Leu_TAG",
    "mttRNALys_TTT": "mt-Lys_TTT",
    "mttRNAMet_CAT": "mt-Met_CAT",
    "mttRNAPhe_GAA": "mt-Phe_GAA",
    "mttRNAPro_TGG": "mt-Pro_TGG",
    "mttRNASer_GCT": "mt-Ser_GCT",
    "mttRNASer_TGA": "mt-Ser_TGA",
    "mttRNAThr_TGT": "mt-Thr_TGT",
    "mttRNATrp_TCA": "mt-Trp_TCA",
    "mttRNATyr_GTA": "mt-Tyr_GTA",
    "mttRNAVal_TAC": "mt-Val_TAC",
})

print(merged_counts_and_gtf_df)

#%%
pd.set_option('display.float_format', '{:.6f}'.format)
source_table = merged_counts_and_gtf_df.groupby("gene_type").sum().drop(columns=["gene_id", "gene_name", "key_0"]).sort_values(by="combined", ascending=False)
source_table = (source_table / 1000000) * 100
print(source_table)
source_table.to_csv(base + "/realigned_bams/gene_type_composition_tRNA_combined_only.tsv", sep="\t",index=True, header=True)


#%%
print(merged_counts_and_gtf_df["gene_type"])
#%%
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 11,          # base
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "legend.title_fontsize": 12,
})

# %%
print(set(merged_counts_and_gtf_df["gene_type"]))
#%%
sample_cols = ["combined"]


comp = merged_counts_and_gtf_df.groupby("gene_type")[sample_cols].sum()   # absolute cpm
rel = comp.div(comp.sum(axis=0), axis=1)                # relative widths

top_n = 8
order = rel.sum(axis=1).sort_values(ascending=False).index
keep = order[:top_n]


def collapse(df):
    out = df.loc[keep].copy()
    out.loc["other"] = df.drop(keep).sum(axis=0)
    return out

rel_plot = collapse(rel)      # controls segment width (0–1)
abs_plot = collapse(comp)     # absolute cpm for labels

fig,ax = plt.subplots()
colors = [gene_type_colors.get(b, "#CCCCCC") for b in rel_plot.index]
rel_plot.T.plot(
    kind="barh", stacked=True, figsize=(16,4),
    color=colors, width=0.8, ax=ax
)

ax.set_xlabel("Relative composition (CPM / 1M)")
ax.set_ylabel("Sample")
ax.set_xlim(0, 1)
ax.legend(title="Biotype", bbox_to_anchor=(1.01, 1), loc="upper left")


plt.tight_layout()
plt.show()
fig.savefig(base + "/realigned_bams/Figure1E_tRNA_biotypes.pdf", format="pdf")

# %%
sample_cols = [
"HRP_A_004_1", 
"HRP_A_006_1",
"HRP_A_012_1",
"HRP_A_020_1_1",
"HRP_A_020_1_2",
"HRP_A_020_1_3",
"HRP_A_020_2_1",
"HRP_A_020_2_2",
"HRP_A_020_2_3",
"HRP_A_020_3_1",
"HRP_A_020_3_2",
"HRP_A_020_3_3",
"HRP_A_021_1_1",
"HRP_A_021_1_2",
"HRP_A_021_1_3"
][::-1]

comp = merged_counts_and_gtf_df.groupby("gene_type")[sample_cols].sum()   # absolute cpm
rel = comp.div(comp.sum(axis=0), axis=1)                # relative widths

# keep gene_types exceeding 2% relative share in at least one sample
keep = rel.index[(rel > 0.01).any(axis=1)]

# order by total signal (so "other" lands last and bars stack consistently)
keep = rel.loc[keep].sum(axis=1).sort_values(ascending=False).index

def collapse(df):
    out = df.loc[keep].copy()
    out.loc["other"] = df.drop(keep).sum(axis=0)
    return out


rel_plot = collapse(rel)      # controls segment width (0–1)
abs_plot = collapse(comp)     # absolute cpm for labels


fig,ax = plt.subplots()
colors = [gene_type_colors.get(b, "#CCCCCC") for b in rel_plot.index]

rel_plot.T.plot(
    kind="barh", stacked=True, figsize=(16,8),
    color=colors, width=0.8, ax=ax
)


ax.set_xlabel("Relative composition (CPM / 1M)")
ax.set_ylabel("tRNA")
ax.set_yticklabels("")
ax.set_xlim(0, 1)
ax.legend(title="Biotype", bbox_to_anchor=(1.01, 1), loc="upper left")

plt.tight_layout()
plt.show()

fig.savefig(base + "/realigned_bams/Figure1_tRNA_biotype_CPMS.pdf", format="pdf")
# %%
distribution_tRNA_df = merged_counts_and_gtf_df[merged_counts_and_gtf_df["gene_type"] == "tRNA"]

print(distribution_tRNA_df.columns)
distribution_tRNA_df['key_0'] = distribution_tRNA_df['key_0'].replace({
    # rRNA
    "hs_rRNA_5S": "5S",
    "hs_rRNA_5.8S": "5.8S",
    "hs_rRNA_18S": "18S",
    "hs_rRNA_28S": "28S",
    "hs_mt-rRNA_12s": "mt-12S",
    "hs_mt-rRNA_16s": "mt-16S",
    #hs_ cytoplasmic tRNA
    "hs_tRNAAla_CGC": "Ala_CGC",
    "hs_tRNAAla_TGC": "Ala_TGC",
    "hs_tRNAAla_AGC": "Ala_AGC",
    "hs_tRNAArg_ACG": "Arg_ACG",
    "hs_tRNAArg_TCG": "Arg_TCG",
    "hs_tRNAArg_CCG2": "Arg_CCG2",
    "hs_tRNAArg_TCT2": "Arg_TCT2",
    "hs_tRNAArg_CCT": "Arg_CCT",
    "hs_tRNAArg_CCG1": "Arg_CCG1",
    "hs_tRNAArg_TCT1": "Arg_TCT1",
    "hs_tRNAAsn_GTT": "Asn_GTT",
    "hs_tRNAAsp_GTC": "Asp_GTC",
    "hs_tRNACys_GCA": "Cys_GCA",
    "hs_tRNAGln_CTG_TTG": "Gln_CTG_TTG",
    "hs_tRNAGlu_CTC": "Glu_CTC",
    "hs_tRNAGlu_TTC": "Glu_TTC",
    "hs_tRNAGly_GCC": "Gly_GCC",
    "hs_tRNAGly_TCC": "Gly_TCC",
    "hs_tRNAGly_CCC": "Gly_CCC",
    "hs_tRNAHis_GTG": "His_GTG",
    "hs_tRNAIle_AAT": "Ile_AAT",
    "hs_tRNAIle_TAT": "Ile_TAT",
    "hs_tRNAIle_GAT": "Ile_GAT",
    "hs_tRNALeu_AAG": "Leu_AAG",
    "hs_tRNALeu_CAG": "Leu_CAG",
    "hs_tRNALeu_CAA": "Leu_CAA",
    "hs_tRNALeu_TAA": "Leu_TAA",
    "hs_tRNALys_CTT": "Lys_CTT",
    "hs_tRNALys_TTT": "Lys_TTT",
    "hs_tRNAMet_CAT": "Met_CAT",
    "hs_tRNAPhe_GAA": "Phe_GAA",
    "hs_tRNAPro_AGG_CGG_TGG": "Pro_AGG_CGG_TGG",
    "hs_tRNASeC_TCA": "SeC_TCA",
    "hs_tRNASer_AGA": "Ser_AGA",
    "hs_tRNASer_CGA": "Ser_CGA",
    "hs_tRNASer_GCT": "Ser_GCT",
    "hs_tRNAThr_AGT": "Thr_AGT",
    "hs_tRNAThr_CGT1": "Thr_CGT1",
    "hs_tRNAThr_TGT": "Thr_TGT",
    "hs_tRNAThr_CGT2": "Thr_CGT2",
    "hs_tRNATrp_CCA": "Trp_CCA",
    "hs_tRNATyr_GTA2": "Tyr_GTA2",
    "hs_tRNATyr_ATA": "Tyr_ATA",
    "hs_tRNATyr_GTA1": "Tyr_GTA1",
    "hs_tRNAVal_AAC_CAC": "Val_AAC_CAC",
    "hs_tRNAVal_TAC": "Val_TAC",
    #hs_ mitochondrial tRNA
    "hs_mttRNAAla_TGC": "mt-Ala_TGC",
    "hs_mttRNAArg_TCG": "mt-Arg_TCG",
    "hs_mttRNAAsn_GTT": "mt-Asn_GTT",
    "hs_mttRNAAsp_GTC": "mt-Asp_GTC",
    "hs_mttRNACys_GCA": "mt-Cys_GCA",
    "hs_mttRNAGln_TTG": "mt-Gln_TTG",
    "hs_mttRNAGlu_TTC": "mt-Glu_TTC",
    "hs_mttRNAGly_TCC": "mt-Gly_TCC",
    "hs_mttRNAHis_GTG": "mt-His_GTG",
    "hs_mttRNAIle_GAT": "mt-Ile_GAT",
    "hs_mttRNALeu_TAA": "mt-Leu_TAA",
    "hs_mttRNALeu_TAG": "mt-Leu_TAG",
    "hs_mttRNALys_TTT": "mt-Lys_TTT",
    "hs_mttRNAMet_CAT": "mt-Met_CAT",
    "hs_mttRNAPhe_GAA": "mt-Phe_GAA",
    "hs_mttRNAPro_TGG": "mt-Pro_TGG",
    "hs_mttRNASer_GCT": "mt-Ser_GCT",
    "hs_mttRNASer_TGA": "mt-Ser_TGA",
    "hs_mttRNAThr_TGT": "mt-Thr_TGT",
    "hs_mttRNATrp_TCA": "mt-Trp_TCA",
    "hs_mttRNATyr_GTA": "mt-Tyr_GTA",
    "hs_mttRNAVal_TAC": "mt-Val_TAC",
})

transcript_colors = {
    # rRNA — ocean blues
    "5S":      "#1F4E79",
    "5.8S":    "#2E6CA4",
    "18S":     "#3D8BC0",
    "28S":     "#5BA3D0",
    "mt-12S":  "#4A6670",   # mt — slate
    "mt-16S":  "#637E88",

    # cytoplasmic tRNA — green→earth amino-acid spectrum
    "Ala_CGC":         "#2E5E3A",
    "Ala_TGC":         "#3A6E45",
    "Ala_AGC":         "#477E50",
    "Arg_ACG":         "#548D5B",
    "Arg_TCG":         "#629C66",
    "Arg_CCG2":        "#5E9B5E",
    "Arg_TCT2":        "#6FA86B",
    "Arg_CCT":         "#7DB073",
    "Arg_CCG1":        "#8AB87B",
    "Arg_TCT1":        "#97C084",
    "Asn_GTT":         "#6B9362",
    "Asp_GTC":         "#759C5A",
    "Cys_GCA":         "#83A551",
    "Gln_CTG_TTG":     "#90AD49",
    "Glu_CTC":         "#9DB541",
    "Glu_TTC":         "#A8BC45",
    "Gly_GCC":         "#B0C04A",
    "Gly_TCC":         "#A3A847",
    "Gly_CCC":         "#959E42",
    "His_GTG":         "#AEB158",
    "Ile_AAT":         "#B8B05A",
    "Ile_TAT":         "#C2AF55",
    "Ile_GAT":         "#CBAE4E",
    "Leu_AAG":         "#D4AD47",
    "Leu_CAG":         "#DAA520",
    "Leu_CAA":         "#D6A032",
    "Leu_TAA":         "#CC9A3A",
    "Lys_CTT":         "#C29440",
    "Lys_TTT":         "#B98C3E",
    "Met_CAT":         "#B5651D",
    "Phe_GAA":         "#AD7A35",
    "Pro_AGG_CGG_TGG": "#A86E2E",
    "SeC_TCA":         "#9C6328",
    "Ser_AGA":         "#A0522D",
    "Ser_CGA":         "#985534",
    "Ser_GCT":         "#8F5A3C",
    "Thr_AGT":         "#8B4F2E",
    "Thr_CGT1":        "#7E5536",
    "Thr_TGT":         "#76502F",
    "Thr_CGT2":        "#6E4A2B",
    "Trp_CCA":         "#7A5A3A",
    "Tyr_GTA2":        "#6B4F38",
    "Tyr_ATA":         "#5F4A33",
    "Tyr_GTA1":        "#54442F",
    "Val_AAC_CAC":     "#4E5340",
    "Val_TAC":         "#5A6049",

    # mitochondrial tRNA — terracotta / rust
    "mt-Ala_TGC":  "#C1693C",
    "mt-Arg_TCG":  "#CD7548",
    "mt-Asn_GTT":  "#D98153",
    "mt-Asp_GTC":  "#B85A30",
    "mt-Cys_GCA":  "#C26538",
    "mt-Gln_TTG":  "#CE7040",
    "mt-Glu_TTC":  "#A84E2A",
    "mt-Gly_TCC":  "#B45932",
    "mt-His_GTG":  "#C0643A",
    "mt-Ile_GAT":  "#9C4A2A",
    "mt-Leu_TAA":  "#A85432",
    "mt-Leu_TAG":  "#B45F3A",
    "mt-Lys_TTT":  "#92442A",
    "mt-Met_CAT":  "#9E4E32",
    "mt-Phe_GAA":  "#AA583A",
    "mt-Pro_TGG":  "#883E28",
    "mt-Ser_GCT":  "#944830",
    "mt-Ser_TGA":  "#A05238",
    "mt-Thr_TGT":  "#7E3826",
    "mt-Trp_TCA":  "#8A422E",
    "mt-Tyr_GTA":  "#964C36",
    "mt-Val_TAC":  "#723424",
}

comp = merged_counts_and_gtf_df.groupby("gene_type")[sample_cols].sum()   # absolute cpm
rel = comp.div(comp.sum(axis=0), axis=1)                # relative widths

# keep gene_types exceeding 2% relative share in at least one sample
keep = rel.index[(rel > 0.01).any(axis=1)]

# order by total signal (so "other" lands last and bars stack consistently)
keep = rel.loc[keep].sum(axis=1).sort_values(ascending=False).index

def collapse(df):
    out = df.loc[keep].copy()
    out.loc["other"] = df.drop(keep).sum(axis=0)
    return out


fig,ax = plt.subplots()
colors = [transcript_colors.get(b, "#CCCCCC") for b in rel_plot.index]

rel_plot.T.plot(
    kind="barh", stacked=True, figsize=(16,8),
    color=colors, width=0.8, ax=ax
)

ax.set_xlabel("Relative composition (CPM / 1M)")
ax.set_ylabel("Sample")
ax.set_xlim(0, 1)
ax.legend(title="Biotype", bbox_to_anchor=(1.01, 1), loc="upper left")

fig.savefig(base + "/realigned_bams/Figure1_tRNA_anticodon_aminoacid_distribution_CPMS.pdf", format="pdf")


# %%
