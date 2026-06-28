process DIVERSITY {
    tag    "DiversityAnalysis"
    label  'diversity'
    publishDir "${params.outdir}/diversity", mode: 'copy', pattern: "*.{html,tsv,png}"

    input:
    path(bracken_reports)

    output:
    path("alpha_diversity.tsv"),    emit: alpha_table
    path("beta_diversity.tsv"),     emit: beta_table
    path("alpha_*.html"),           emit: alpha_plots
    path("beta_*.html"),            emit: beta_plots
    path("ordination_*.html"),      emit: ordination_plots

    script:
    def reports_json = groovy.json.JsonOutput.toJson(bracken_reports.collect { it.toString() })
    """
    # Write file list for Python
    cat > .reports.json << 'JSONEOF'
${reports_json}
JSONEOF

    python3 <<CODE
import os, json, warnings
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE
from sklearn.preprocessing import StandardScaler
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

warnings.filterwarnings('ignore')

with open('.reports.json') as f:
    reports = json.load(f)
if not reports:
    print("No Bracken reports found. Creating empty outputs.")
    for f in ["alpha_diversity.tsv", "beta_diversity.tsv",
              "alpha_boxplots.html", "alpha_violin.html",
              "beta_heatmap.html", "ordination_pca.html",
              "ordination_pcoa.html", "ordination_nmds.html",
              "ordination_umap.html", "ordination_tsne.html"]:
        with open(f, 'w') as fh:
            fh.write("")
    exit(0)

# Parse all Bracken species reports into an abundance matrix
samples = []
abundances = []
for r in sorted(reports):
    sid = os.path.basename(r).replace('.bracken.S.report', '').replace('.bracken', '')
    df = pd.read_csv(r, sep='\\t')
    df.columns = [c.strip() for c in df.columns]
    df = df.set_index('name')
    samples.append(sid)
    abundances.append(df['new_est_reads'])

abund = pd.concat(abundances, axis=1, join='outer').fillna(0)
abund.columns = samples
abund = abund.loc[(abund.sum(axis=1) > 0)]
n_species = abund.shape[0]

# Normalize to relative abundance (%)
rel_abund = abund.div(abund.sum(axis=0), axis=1) * 100

# ── 1. Alpha Diversity ──────────────────────────────────────────────────
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
    return float(obs)

alpha_df = pd.DataFrame({
    'Sample': abund.columns,
    'Shannon': abund.apply(shannon),
    'Simpson': abund.apply(simpson),
    'Chao1': abund.apply(chao1),
    'Observed': (abund > 0).sum(),
})
alpha_df.to_csv("alpha_diversity.tsv", sep='\\t', index=False)

# ── Alpha Boxplots ──────────────────────────────────────────────────────
metrics = ['Shannon', 'Simpson', 'Chao1', 'Observed']
fig_box = make_subplots(rows=2, cols=2,
    subplot_titles=metrics,
    vertical_spacing=0.15, horizontal_spacing=0.12)

for i, m in enumerate(metrics):
    row, col = (i // 2) + 1, (i % 2) + 1
    fig_box.add_trace(
        go.Box(y=alpha_df[m], name=m,
               marker_color='#2E86AB', boxmean='sd',
               hovertemplate='%{y:.3f}<extra></extra>'),
        row=row, col=col)

fig_box.update_layout(title='Alpha Diversity Metrics — Boxplots',
                      height=600, showlegend=False,
                      template='plotly_white')
fig_box.write_html("alpha_boxplots.html")

# ── Alpha Violin Plots ──────────────────────────────────────────────────
fig_violin = make_subplots(rows=2, cols=2,
    subplot_titles=metrics,
    vertical_spacing=0.15, horizontal_spacing=0.12)

for i, m in enumerate(metrics):
    row, col = (i // 2) + 1, (i % 2) + 1
    fig_violin.add_trace(
        go.Violin(y=alpha_df[m], name=m,
                  box_visible=True, meanline_visible=True,
                  line_color='#2E86AB', fillcolor='#A2D5F2',
                  opacity=0.7,
                  hovertemplate='%{y:.3f}<extra></extra>'),
        row=row, col=col)

fig_violin.update_layout(title='Alpha Diversity Metrics — Violin Plots',
                         height=600, showlegend=False,
                         template='plotly_white')
fig_violin.write_html("alpha_violin.html")

# ── 2. Beta Diversity ───────────────────────────────────────────────────
# Bray-Curtis
bc = pdist(rel_abund.T, metric='braycurtis')
bc_mat = squareform(bc)
labels = list(rel_abund.columns)

# Jaccard
jacc = pdist(rel_abund.T > 0, metric='jaccard')
jacc_mat = squareform(jacc)

beta_dict = {'Bray-Curtis': bc_mat, 'Jaccard': jacc_mat}

# Beta diversity heatmaps
for metric_name, mat in beta_dict.items():
    fig = go.Figure(data=go.Heatmap(
        z=mat, x=labels, y=labels,
        colorscale='YlOrRd', zmin=0, zmax=1,
        text=np.round(mat, 3),
        texttemplate='%{text}',
        textfont=dict(size=8),
        hovertemplate='%{x} vs %{y}<br>%{metric_name}: %{z:.3f}<extra></extra>'))
    fig.update_layout(
        title=f'Beta Diversity — {metric_name} Dissimilarity',
        xaxis_title='Sample', yaxis_title='Sample',
        height=500, width=550, template='plotly_white')
    fname = f"beta_{metric_name.lower().replace(' ', '_')}.html"
    fig.write_html(fname)

# Combined beta table
beta_out = pd.DataFrame({
    'Sample_1': np.repeat(labels, len(labels)),
    'Sample_2': np.tile(labels, len(labels)),
    'Bray_Curtis': bc_mat.flatten(),
    'Jaccard': jacc_mat.flatten(),
})
beta_out.to_csv("beta_diversity.tsv", sep='\\t', index=False)

# Hierarchical clustering dendrogram based on Bray-Curtis
Z = linkage(bc_mat, method='average')
fig_dend = go.Figure()
fig_dend.add_trace(go.Scatter(
    x=[], y=[], mode='markers',
    marker=dict(size=0)))
fig_dend.update_layout(
    title='Hierarchical Clustering (Bray-Curtis, UPGMA)',
    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    yaxis=dict(title='Distance'),
    height=400, template='plotly_white')
# Add dendrogram using shapes
from plotly.figure_factory import create_dendrogram
dendro = create_dendrogram(rel_abund.T, orientation='bottom',
                           labels=labels, linkagefun=lambda x: Z)
dendro.update_layout(title='Hierarchical Clustering (Bray-Curtis, UPGMA)',
                     height=400, template='plotly_white')
dendro.write_html("beta_dendrogram.html")

# ── 3. Ordinations ──────────────────────────────────────────────────────
# PCA
scaler = StandardScaler()
X = scaler.fit_transform(rel_abund.T)
pca = PCA(n_components=min(4, X.shape[0], X.shape[1]))
coords = pca.fit_transform(X)
var_exp = pca.explained_variance_ratio_ * 100

fig_pca = go.Figure()
colors = px.colors.qualitative.Set2
for i, sid in enumerate(labels):
    fig_pca.add_trace(go.Scatter(
        x=[coords[i, 0]], y=[coords[i, 1]],
        mode='markers+text',
        text=[sid], textposition='top center',
        marker=dict(size=14, color=colors[i % len(colors)],
                    line=dict(width=1, color='black')),
        name=sid,
        hovertemplate='%{text}<br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>'))
fig_pca.update_layout(
    title=f'PCA — {var_exp[0]:.1f}% / {var_exp[1]:.1f}% variance',
    xaxis_title=f'PC1 ({var_exp[0]:.1f}%)',
    yaxis_title=f'PC2 ({var_exp[1]:.1f}%)',
    height=500, width=650, template='plotly_white',
    hovermode='closest')
fig_pca.write_html("ordination_pca.html")

# PCoA (MDS) on Bray-Curtis
mds = MDS(n_components=2, dissimilarity='precomputed',
          random_state=42, normalized_stress=False)
mds_coords = mds.fit_transform(bc_mat)

fig_pcoa = go.Figure()
for i, sid in enumerate(labels):
    fig_pcoa.add_trace(go.Scatter(
        x=[mds_coords[i, 0]], y=[mds_coords[i, 1]],
        mode='markers+text',
        text=[sid], textposition='top center',
        marker=dict(size=14, color=colors[i % len(colors)],
                    line=dict(width=1, color='black')),
        name=sid,
        hovertemplate='%{text}<br>MDS1: %{x:.2f}<br>MDS2: %{y:.2f}<extra></extra>'))
fig_pcoa.update_layout(
    title='PCoA (MDS on Bray-Curtis)',
    xaxis_title='MDS1', yaxis_title='MDS2',
    height=500, width=650, template='plotly_white',
    hovermode='closest')
fig_pcoa.write_html("ordination_pcoa.html")

# NMDS
from sklearn.manifold import MDS as SkMDS
class NMDS:
    def __init__(self, n_components=2, random_state=42):
        self.mds = SkMDS(n_components=n_components,
                         dissimilarity='precomputed',
                         random_state=random_state,
                         normalized_stress=False)
    def fit_transform(self, D):
        return self.mds.fit_transform(D)

nmds = NMDS(n_components=2)
nmds_coords = nmds.fit_transform(bc_mat)

fig_nmds = go.Figure()
for i, sid in enumerate(labels):
    fig_nmds.add_trace(go.Scatter(
        x=[nmds_coords[i, 0]], y=[nmds_coords[i, 1]],
        mode='markers+text',
        text=[sid], textposition='top center',
        marker=dict(size=14, color=colors[i % len(colors)],
                    line=dict(width=1, color='black')),
        name=sid,
        hovertemplate='%{text}<br>NMDS1: %{x:.2f}<br>NMDS2: %{y:.2f}<extra></extra>'))
fig_nmds.update_layout(
    title='NMDS (Non-metric Multidimensional Scaling, Bray-Curtis)',
    xaxis_title='NMDS1', yaxis_title='NMDS2',
    height=500, width=650, template='plotly_white',
    hovermode='closest')
fig_nmds.write_html("ordination_nmds.html")

# UMAP
try:
    import umap
    reducer = umap.UMAP(n_components=2, random_state=42, metric='cosine')
    umap_coords = reducer.fit_transform(rel_abund.T.values)

    fig_umap = go.Figure()
    for i, sid in enumerate(labels):
        fig_umap.add_trace(go.Scatter(
            x=[umap_coords[i, 0]], y=[umap_coords[i, 1]],
            mode='markers+text',
            text=[sid], textposition='top center',
            marker=dict(size=14, color=colors[i % len(colors)],
                        line=dict(width=1, color='black')),
            name=sid,
            hovertemplate='%{text}<br>UMAP1: %{x:.2f}<br>UMAP2: %{y:.2f}<extra></extra>'))
    fig_umap.update_layout(
        title='UMAP (cosine distance on relative abundance)',
        xaxis_title='UMAP1', yaxis_title='UMAP2',
        height=500, width=650, template='plotly_white',
        hovermode='closest')
    fig_umap.write_html("ordination_umap.html")
except ImportError:
    with open("ordination_umap.html", 'w') as f:
        f.write('<html><body><h3>UMAP not available (install umap-learn)</h3></body></html>')

# t-SNE
tsne = TSNE(n_components=2, random_state=42, metric='euclidean',
            perplexity=min(30, len(labels) - 1) if len(labels) > 1 else 1)
tsne_coords = tsne.fit_transform(rel_abund.T.values)

fig_tsne = go.Figure()
for i, sid in enumerate(labels):
    fig_tsne.add_trace(go.Scatter(
        x=[tsne_coords[i, 0]], y=[tsne_coords[i, 1]],
        mode='markers+text',
        text=[sid], textposition='top center',
        marker=dict(size=14, color=colors[i % len(colors)],
                    line=dict(width=1, color='black')),
        name=sid,
        hovertemplate='%{text}<br>t-SNE1: %{x:.2f}<br>t-SNE2: %{y:.2f}<extra></extra>'))
fig_tsne.update_layout(
    title='t-SNE (Euclidean on relative abundance)',
    xaxis_title='t-SNE1', yaxis_title='t-SNE2',
    height=500, width=650, template='plotly_white',
    hovermode='closest')
fig_tsne.write_html("ordination_tsne.html")

print(f"Diversity analysis complete: {len(labels)} samples, {n_species} species")
CODE
    """

    stub:
    """
    touch alpha_diversity.tsv beta_diversity.tsv \
          alpha_boxplots.html alpha_violin.html \
          beta_bray_curtis.html beta_jaccard.html beta_dendrogram.html \
          ordination_pca.html ordination_pcoa.html ordination_nmds.html \
          ordination_umap.html ordination_tsne.html
    """
}
