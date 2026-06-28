#!/usr/bin/env python3
import sys
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

sample = sys.argv[1]
bracken_p = sys.argv[2]
bracken_g = sys.argv[3]
bracken_s = sys.argv[4]

def read_bracken(path):
    df = pd.read_csv(path, sep='\t')
    df.columns = [c.strip() for c in df.columns]
    return df

try:
    df_p = read_bracken(bracken_p)
    df_g = read_bracken(bracken_g)
    df_s = read_bracken(bracken_s)
except Exception as e:
    print(f"No Bracken data for {sample}, creating empty plots")
    open(f"{sample}.sankey.png", 'w').write('')
    open(f"{sample}.summary.png", 'w').write('')
    exit(0)

top_phyla = df_p.nlargest(10, 'new_est_reads')
top_phyla_names = set(top_phyla['name'])

df_g_top = df_g[df_g['name'].isin(top_phyla_names)]
if len(df_g_top) < 5:
    df_g_top = df_g.nlargest(20, 'new_est_reads')

top_genera = df_g_top.nlargest(15, 'new_est_reads')

sources = []
targets = []
values = []
node_labels = []
node_map = {}

for _, row in top_phyla.iterrows():
    name = row['name']
    if name not in node_map:
        node_map[name] = len(node_labels)
        node_labels.append(name)

for _, row in top_genera.iterrows():
    name = row['name']
    if name not in node_map:
        node_map[name] = len(node_labels)
        node_labels.append(name)

genus_phyla = {}
for _, prow in top_phyla.iterrows():
    pname = prow['name']
    reads = prow['new_est_reads']

for _, grow in top_genera.iterrows():
    gname = grow['name']
    greads = grow['new_est_reads']
    best_phy = None
    best_score = 0
    for pname in top_phyla_names:
        if pname.lower() in gname.lower() or gname.lower() in pname.lower():
            score = 1
            if score > best_score:
                best_score = score
                best_phy = pname
    if best_phy is None:
        best_phy = top_phyla.iloc[0]['name']
    if greads > 0:
        sources.append(node_map[best_phy])
        targets.append(node_map[gname])
        values.append(greads)

if not values:
    sources.append(0)
    targets.append(0)
    values.append(1)
    node_labels = ['No data']

fig_sankey = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color='black', width=0.5),
        label=node_labels,
    ),
    link=dict(
        source=sources,
        target=targets,
        value=values,
    ))])

fig_sankey.update_layout(
    title=f'{sample} — Phylum to Genus Flow (Sankey)',
    font_size=12,
    height=600,
)
fig_sankey.write_image(f"{sample}.sankey.png", format='png', width=900, height=600)

fig = make_subplots(
    rows=2, cols=2,
    subplot_titles=('Top Phyla', 'Top Genera', 'Top Species', 'Classification'),
    vertical_spacing=0.15,
)

top10_p = df_p.nlargest(10, 'new_est_reads')
fig.add_trace(
    go.Bar(y=top10_p['name'][::-1], x=top10_p['new_est_reads'][::-1],
           orientation='h', marker_color='#2E86AB', name='Phyla'),
    row=1, col=1
)

top10_g = df_g.nlargest(10, 'new_est_reads')
fig.add_trace(
    go.Bar(y=top10_g['name'][::-1], x=top10_g['new_est_reads'][::-1],
           orientation='h', marker_color='#A23B72', name='Genera'),
    row=1, col=2
)

top10_s = df_s.nlargest(10, 'new_est_reads')
fig.add_trace(
    go.Bar(y=top10_s['name'][::-1], x=top10_s['new_est_reads'][::-1],
           orientation='h', marker_color='#F18F01', name='Species'),
    row=2, col=1
)

fig.update_layout(
    title=f'{sample} — Taxonomic Summary',
    height=800,
    showlegend=False,
)
fig.write_image(f"{sample}.summary.png", format='png', width=900, height=800)
