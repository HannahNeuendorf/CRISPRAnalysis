# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 14:16:19 2022

@author: hannahN
"""
## This script must be run in Python 3.9  ##

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

###Plot by hit/non-hit 

# Load data
volcano_df = pd.read_csv("MM576_volcano_AdhCTRLvNonAdh_R3.csv", sep=",", index_col=0)
#print(interested_genes)

# Apply negative log10 transform
volcano_df["-p-value"] = -np.log10(
    volcano_df["p-value"]
)

# Filter out pseudo genes
volcano_df = volcano_df[
    volcano_df.index.isin([i for i in volcano_df.index if i[:6] != "pseudo"])
]

# p-value thresholding
pvalue_thresthold = -np.log10(0.05)
volcano_df["Hit/Non-hit"] = "Non-hit"
volcano_df.loc[
    volcano_df["-p-value"] >= pvalue_thresthold,
    ["Hit/Non-hit"],
] = "Hit"

# Switch non/non-hit labels for genes of interest to "Genes of interest"
interested_genes = ["DYNLRB1", "OR51B4", "EPS8L1", "PRSS50", "DUSP14"]
volcano_df.loc[
    volcano_df.index.isin(interested_genes), ["Hit/Non-hit"]
] = "Genes of interest"

color_map_d = {"Non-hit": "lightgray", "Hit": "red", "Genes of interest": "green"}


#Get list of hit genes 
hit_genes = volcano_df[volcano_df["Hit/Non-hit"] == "Hit"].index.tolist()

# Get list of genes as a separate coluns
volcano_df["genes"] = volcano_df.index

fig = px.scatter(
    volcano_df,
    x="log2FC",
    y="-p-value",
    title="MM576 AdhCTRL v NonAdh Rep3",
    color="Hit/Non-hit",
    hover_name="genes",
    hover_data=["log2FC", "p-value"],
       color_discrete_map=color_map_d,
    category_orders={"Hit/Non-hit": ["Hit", "Non-hit"]},
)


# Update scatter plot properties
fig.update_traces(opacity=1, marker=dict(size=2))

# Make Genes of interest disappear in legend
fig.update_traces(selector=dict(name="Genes of interest"), showlegend=False)

# Annotate genes of interest
for gene in interested_genes:
    x_coor = volcano_df.loc[gene, "log2FC"]
    y_coor = volcano_df.loc[gene, "-p-value"]
    fig.add_annotation(
        x=x_coor,
        y=y_coor,
        text=gene,
        textangle=0,
        font_size=15,
        showarrow=False,
        yshift=0,
        xshift=40,
    )

# Update scatter plot properties
fig.update_traces(opacity=1, marker=dict(size=3))

fig.update_traces(markermarker=dict(color="green"), selector=dict(type="Highlighed"))

# Add significant logfold of p-value line
# Comment this out if you don't want lines
fig.add_hline(
    y=pvalue_thresthold,
    line_dash="dash",
    line_color="gray",
    line_width=0.5,
    # annotation_text="Jan 1, 2018 baseline",
    # annotation_position="bottom right",
)

# Adjust axes
fig.update_yaxes(
    title="-log10 p-value",
    title_font_size=12,
    title_standoff=10,
    linecolor="black",
    linewidth=1,
    ticks="outside",
    ticklen=2,
    tickwidth=1,
    dtick=0.2,
    range=[0, 3.1],
    tickfont_size=8,
    # showgrid=True,
    # gridcolor="lightgray"
)
fig.update_xaxes(
    title="log2FC",
    title_font_size=12,
    title_standoff=5,
    linecolor="black",
    linewidth=1,
    ticks="outside",
    ticklen=2,
    tickwidth=0.5,
    tickfont_size=8,
    dtick=2.5,
    range=[-15, 21],
    # showgrid=True,
    # gridcolor="lightgray"
)

# Update layout
fig["layout"].update(
    font=dict(
        size=10,
    ),
    font_family="Times New Roman",
    plot_bgcolor="rgba(0,0,0,0)",
    legend=dict(title="Gene types", title_font_size=14, font_size=12),
    showlegend=True,
    newshape=dict(opacity=1),
    margin=dict(t=20, l=20, r=0, b=0),  # Tight margin
)

fig.write_image("MM576_AdhCTRLvNonAdh_Rep3_volcano.pdf", scale=5, width=600, height=400)


# Save into HTML
fig.write_html(
    ("MM576_AdhCTRLvNonAdh_Rep3_volcano.html"), auto_open=True,
)

