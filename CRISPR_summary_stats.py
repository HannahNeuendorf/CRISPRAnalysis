# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 08:21:57 2022

@author: hannahN
"""
import numpy as np
import pandas as pd

import plotly
import plotly.express as px

from pathlib import Path
##note before loading in a csv file, need to combine the top two row labels, add sgRNA in A1, save as .csv and remove the extra space between sgRNA and the comma
# Set prefix
#prefix = "/SequencingData/20220815_NovaSeq_370_576_CRISPRa/process_experiments/MM576_output"
# Read data files
#df = pd.read_csv(
#    Path(prefix).joinpath("CRISPRa_MM370_mergedcountstable.csv"),
#    index_col=0,
#    sep=","
#)

df = pd.read_csv(
    ("CRISPRa_MM576_mergedcountstable_NonAdh.csv"),
    index_col=0,
    sep=","
    )


# Column names are samples
samples = df.columns.tolist()

# Reset index to turn sgRNA into columN
df.reset_index(inplace=True)

# Extract gene symbols
df["gene_symbols"] = [i.split("_")[0] for i in df["sgRNA"]]

# Move gene symbols to first columns
gene_sym_col = df.pop("gene_symbols")
df.insert(0, "gene_symbols", gene_sym_col)

# Filter out controls
df = df[
    df["gene_symbols"]!="non-targeting"
]

##Tally RNA counts by sgRNA
# Calculate total counts by samples and genes
total_counts_by_sample_by_genes_df = df.groupby(["gene_symbols"])[samples].sum()

# Calculate number of identified genes
identified_genes = (total_counts_by_sample_by_genes_df.sum(axis=1)!=0).sum()

# Calculate number of missing genes
total_genes = 18915
missing_genes = total_genes - identified_genes

# Make a DataFrame hosting number of missing and identified genes
plot_rna_df = pd.DataFrame(
    columns=["Identified/Absent", "Gene_counts"],
    data=[
        ["Identified", identified_genes],
        ["Absent", missing_genes],
    ]
)
plot_rna_df


# Plot pie chart beautifully
fig = px.pie(
    plot_rna_df, 
    values='Gene_counts', 
    names='Identified/Absent',
    #title="Number of identified & missing genes",
    labels={"Identified/Absent": "Gene_counts"},
    hole=0.4
)

# Update traces and layout
fig.update_traces(textposition='outside', textinfo='percent+label')

fig.update_layout(
    margin=dict(t=0, l=0, r=0, b=0),
    font=dict(size=10),
    font_family="Times New Roman",
    showlegend=True
)

# Save into HTML
fig.write_html(
    ("MM576_NonAdh_genes_covered.html"), auto_open=False
)

# Save into image
fig.write_image(
    ("MM576_NonAdh_genes_covered.pdf"), scale=5, width=350, height=300
)




##Tally RNA counts by gene symbols 
df.groupby(["sgRNA"])[samples].sum().sum(axis=1).to_frame().rename(
    columns={0: "total_counts"}
)

# Calculate total counts by samples and guides
total_counts_by_sample_by_guides_df = df.groupby(["sgRNA"])[samples].sum()

# Calculate number of identified guides
identified_guides = (total_counts_by_sample_by_guides_df.sum(axis=1)!=0).sum()

# Calculate number of missing guides
total_guides = 104540-1895
missing_guides = total_guides - identified_guides

# Make a DataFrame hosting number of missing and identified gens
plot_guides_df = pd.DataFrame(
    columns=["Identified/Absent", "Guides_counts"],
    data=[
        ["Identified", identified_guides],
        ["Absent", missing_guides],
    ]
)
plot_guides_df
# Plot pie chart beautifully
fig = px.pie(
    plot_guides_df, 
    values='Guides_counts', 
    names='Identified/Absent',
    #title="Number of identified & missing sgRNAs",
    labels={"Identified/Absent": "Guides_counts"},
    hole=0.4
)

# Update traces and layout
fig.update_traces(textposition='outside', textinfo='percent+label')

fig.update_layout(
    margin=dict(t=0, l=0, r=0, b=0),
    font=dict(size=10),
    font_family="Times New Roman",
    showlegend=False
)

# Save into HTML
fig.write_html(
    ("MM576_NonAdh_sgRNA_covered.html"), auto_open=False
)

# Save into image - can change svg to any format you like
fig.write_image(
    ("MM576_NonAdh_sgRNA_covered.pdf"), scale=5, width=350, height=300,
    )