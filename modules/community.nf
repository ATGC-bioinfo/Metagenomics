process COMMUNITY {
    tag    "CommunityAnalysis"
    label  'community'
    publishDir "${params.outdir}/community", mode: 'copy', pattern: "*.html"

    input:
    path(bracken_reports)

    output:
    path("community_alpha.html"), emit: alpha
    path("community_beta.html"),  emit: beta
    path("community_heatmap.html"), emit: heatmap
    path("community_pca.html"),   emit: pca

    script:
    """
    python3 <<CODE
import os, re, glob
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

reports = glob.glob("${bracken_reports}")
samples = []
abundances = []

for r in sorted(reports):
    sid = os.path.basename(r).split('.')[0]
    df = pd.read_csv(r, sep='\\t')
    df.columns = [c.strip() for c in df.columns]
    df = df.set_index('name')
    samples.append(sid)
    abundances.append(df['new_est_reads'])

abund = pd.concat(abundances, axis=1, join='outer').fillna(0)
abund.columns = samples
abund = abund.loc[(abund.sum(axis=1) > 0)]

# Normalize to relative abundance
rel_abund = abund.div(abund.sum(axis=0), axis=1) * 100

# ── Alpha Diversity ──────────────────────────────────────────
def shannon(v):
    v = v[v > 0]
    p = v / v.sum()
    return -np.sum(p * np.log2(p))

def simpson(v):
    v = v[v > 0]
    p = v / v.sum()
    return 1 - np.sum(p ** 2)

def chao1(v):
    obs = (v > 0).sum()
    singletons = (v == 1).sum()
    doubletons = (v == 2).sum()
    if doubletons > 0:
        return obs + (singletons ** 2) / (2 * doubletons)
    return obs

alpha_df = pd.DataFrame({
    'Sample': abund.columns,
    'Shannon': abund.apply(shannon),
    'Simpson': abund.apply(simpson),
    'Chao1': abund.apply(chao1),
    'Observed': (abund > 0).sum(),
})

fig_alpha = make_subplots(rows=1, cols=3,
    subplot_titles=('Shannon Index', 'Simpson Index', 'Chao1 Richness'))
for i, metric in enumerate(['Shannon', 'Simpson', 'Chao1'], 1):
    fig_alpha.add_trace(
        go.Bar(x=alpha_df['Sample'], y=alpha_df[metric],
               marker_color='#2E86AB', name=metric),
        row=1, col=i)
fig_alpha.update_layout(title='Alpha Diversity Metrics', height=400)
fig_alpha.write_html("community_alpha.html")

# ── Beta Diversity (Bray-Curtis) ─────────────────────────────
from scipy.spatial.distance import pdist, squareform
bc = pdist(rel_abund.T, metric='braycurtis')
bc_mat = squareform(bc)
labels = rel_abund.columns

fig_beta = go.Figure(data=go.Heatmap(
    z=bc_mat, x=labels, y=labels,
    colorscale='YlOrRd', zmin=0, zmax=1,
    hovertemplate='%{x} vs %{y}<br>Bray-Curtis: %{z:.3f}<extra></extra>'))
fig_beta.update_layout(
    title='Beta Diversity — Bray-Curtis Dissimilarity',
    xaxis_title='Sample', yaxis_title='Sample',
    height=500, width=500)
fig_beta.write_html("community_beta.html")

# ── PCA ──────────────────────────────────────────────────────
scaler = StandardScaler()
X = scaler.fit_transform(rel_abund.T)
pca = PCA(n_components=min(3, len(samples) - 1, X.shape[1]))
coords = pca.fit_transform(X)
var_exp = pca.explained_variance_ratio_ * 100

fig_pca = go.Figure()
fig_pca.add_trace(go.Scatter(
    x=coords[:, 0], y=coords[:, 1],
    mode='markers+text',
    text=labels,
    marker=dict(size=12, color=list(range(len(labels))),
                colorscale='Viridis', showscale=False),
    hovertemplate='%{text}<br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>'))
fig_pca.update_layout(
    title=f'PCA — {var_exp[0]:.1f}% / {var_exp[1]:.1f}% variance',
    xaxis_title=f'PC1 ({var_exp[0]:.1f}%)',
    yaxis_title=f'PC2 ({var_exp[1]:.1f}%)',
    height=500, width=600)
fig_pca.write_html("community_pca.html")

# ── Heatmap (Top 20 species) ─────────────────────────────────
top20 = rel_abund.mean(axis=1).nlargest(20).index
heat_data = rel_abund.loc[top20]

fig_heat = go.Figure(data=go.Heatmap(
    z=heat_data.values,
    x=heat_data.columns,
    y=heat_data.index,
    colorscale='Viridis',
    hovertemplate='%{y}<br>%{x}: %{z:.2f}%<extra></extra>'))
fig_heat.update_layout(
    title='Top 20 Species — Relative Abundance (%)',
    xaxis_title='Sample',
    yaxis_title='Species',
    height=max(400, 20 * 20),
    margin=dict(l=200))
fig_heat.write_html("community_heatmap.html")
CODE
    """

    stub:
    """
    touch community_alpha.html community_beta.html community_heatmap.html community_pca.html
    """
}
