#!/usr/bin/env python3
import os, warnings
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
from plotly.figure_factory import create_dendrogram

warnings.filterwarnings('ignore')

reports = sorted([f for f in os.listdir('.') if f.endswith('.bracken.S.report')])
if not reports:
    print("No Bracken reports found. Creating empty outputs.")
    for f in ["alpha_diversity.tsv", "beta_diversity.tsv",
              "alpha_boxplots.png", "alpha_violin.png",
              "beta_heatmap.png", "ordination_pca.png",
              "ordination_pcoa.png", "ordination_nmds.png",
              "ordination_umap.png", "ordination_tsne.png"]:
        with open(f, 'w') as fh:
            fh.write("")
    exit(0)

samples = []
abundances = []
for r in sorted(reports):
    sid = os.path.basename(r).replace('.bracken.S.report', '').replace('.bracken', '')
    df = pd.read_csv(r, sep='\t')
    df.columns = [c.strip() for c in df.columns]
    df = df.set_index('name')
    samples.append(sid)
    abundances.append(df['new_est_reads'])

abund = pd.concat(abundances, axis=1, join='outer').fillna(0)
abund.columns = samples
abund = abund.loc[(abund.sum(axis=1) > 0)]
n_species = abund.shape[0]

rel_abund = abund.div(abund.sum(axis=0), axis=1) * 100

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
alpha_df.to_csv("alpha_diversity.tsv", sep='\t', index=False)

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
fig_box.write_image("alpha_boxplots.png", format='png', width=800, height=500)

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
fig_violin.write_image("alpha_violin.png", format='png', width=800, height=500)

bc = pdist(rel_abund.T, metric='braycurtis')
bc_mat = squareform(bc)
labels = list(rel_abund.columns)

jacc = pdist(rel_abund.T > 0, metric='jaccard')
jacc_mat = squareform(jacc)

beta_dict = {'Bray-Curtis': bc_mat, 'Jaccard': jacc_mat}

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
    fname = f"beta_{metric_name.lower().replace(' ', '_')}.png"
    fig.write_image(fname, format='png', width=700, height=600)

beta_out = pd.DataFrame({
    'Sample_1': np.repeat(labels, len(labels)),
    'Sample_2': np.tile(labels, len(labels)),
    'Bray_Curtis': bc_mat.flatten(),
    'Jaccard': jacc_mat.flatten(),
})
beta_out.to_csv("beta_diversity.tsv", sep='\t', index=False)

Z = linkage(bc_mat, method='average')
dendro = create_dendrogram(rel_abund.T, orientation='bottom',
                           labels=labels, linkagefun=lambda x: Z)
dendro.update_layout(title='Hierarchical Clustering (Bray-Curtis, UPGMA)',
                     height=400, template='plotly_white')
dendro.write_image("beta_dendrogram.png", format='png', width=700, height=500)

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
fig_pca.write_image("ordination_pca.png", format='png', width=700, height=500)

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
fig_pcoa.write_image("ordination_pcoa.png", format='png', width=700, height=500)

class NMDS:
    def __init__(self, n_components=2, random_state=42):
        self.mds = MDS(n_components=n_components,
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
fig_nmds.write_image("ordination_nmds.png", format='png', width=700, height=500)

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
    fig_umap.write_image("ordination_umap.png", format='png', width=700, height=500)
except ImportError:
    with open("ordination_umap.png", 'w') as f:
        f.write('')

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
fig_tsne.write_image("ordination_tsne.png", format='png', width=700, height=500)

print(f"Diversity analysis complete: {len(labels)} samples, {n_species} species")
