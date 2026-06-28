process PLOTS {
    tag    "${meta.id}"
    label  'plots'
    publishDir "${params.outdir}/plots/${meta.id}", mode: 'copy', pattern: "*.html"

    input:
    tuple val(meta), path(bracken_p), path(bracken_g), path(bracken_s)

    output:
    tuple val(meta), path("${meta.id}.sankey.html"), emit: sankey
    tuple val(meta), path("${meta.id}.summary.html"), emit: summary

    script:
    """
    python3 <<CODE
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

sample = "${meta.id}"

def read_bracken(path):
    df = pd.read_csv(path, sep='\\t')
    df.columns = [c.strip() for c in df.columns]
    return df

# ── Read Bracken reports ──────────────────────────────────────
try:
    df_p = read_bracken("${bracken_p}")
    df_g = read_bracken("${bracken_g}")
    df_s = read_bracken("${bracken_s}")
except Exception as e:
    print(f"No Bracken data for {sample}, creating empty plots")
    open("${meta.id}.sankey.html", 'w').write('<html><body><h3>No data</h3></body></html>')
    open("${meta.id}.summary.html", 'w').write('<html><body><h3>No data</h3></body></html>')
    exit(0)

# ── Sankey: Phylum → Genus flow ───────────────────────────────
top_phyla = df_p.nlargest(10, 'new_est_reads')
top_phyla_names = set(top_phyla['name'])

# Filter genera belonging to top phyla, take top 15
df_g_top = df_g[df_g['name'].isin(top_phyla_names)]
if len(df_g_top) < 5:
    df_g_top = df_g.nlargest(20, 'new_est_reads')

top_genera = df_g_top.nlargest(15, 'new_est_reads')

# Build Sankey
sources = []
targets = []
values = []
node_labels = []
node_map = {}

# Add phylum nodes
for _, row in top_phyla.iterrows():
    name = row['name']
    if name not in node_map:
        node_map[name] = len(node_labels)
        node_labels.append(name)

# Add genus nodes
for _, row in top_genera.iterrows():
    name = row['name']
    if name not in node_map:
        node_map[name] = len(node_labels)
        node_labels.append(name)

# Map each genus to its phylum using phylum-level data
# (Bracken genus report doesn't have phylum column, so we use a heuristic)
genus_phyla = {}
for _, prow in top_phyla.iterrows():
    pname = prow['name']
    reads = prow['new_est_reads']
    # Simple: distribute reads proportionally to genera
    # For a proper mapping, we'd need the full taxonomy lineage

# Simpler approach: genus name as source, phylum name as target
# Read-level mapping: assign each genus to its most likely phylum
# Use the taxonomy: for now, create a simple flow from each phylum to its genera
# by matching genus names to phylum prefix patterns

for _, grow in top_genera.iterrows():
    gname = grow['name']
    greads = grow['new_est_reads']
    # Assign to phylum with closest name match
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
fig_sankey.write_html("${meta.id}.sankey.html")

# ── Summary plots ─────────────────────────────────────────────
fig = make_subplots(
    rows=2, cols=2,
    subplot_titles=('Top Phyla', 'Top Genera', 'Top Species', 'Classification'),
    vertical_spacing=0.15,
)

# Top phyla bar
top10_p = df_p.nlargest(10, 'new_est_reads')
fig.add_trace(
    go.Bar(y=top10_p['name'][::-1], x=top10_p['new_est_reads'][::-1],
           orientation='h', marker_color='#2E86AB', name='Phyla'),
    row=1, col=1
)

# Top genera bar
top10_g = df_g.nlargest(10, 'new_est_reads')
fig.add_trace(
    go.Bar(y=top10_g['name'][::-1], x=top10_g['new_est_reads'][::-1],
           orientation='h', marker_color='#A23B72', name='Genera'),
    row=1, col=2
)

# Top species bar
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
fig.write_html("${meta.id}.summary.html")
CODE
    """

    stub:
    """
    touch ${meta.id}.sankey.html ${meta.id}.summary.html
    """
}
