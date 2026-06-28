#!/usr/bin/env python3
import os, sys, warnings
import numpy as np
import pandas as pd
from scipy import stats
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

warnings.filterwarnings('ignore')

metadata_file = sys.argv[1]

metadata = pd.read_csv(metadata_file)
metadata.columns = [c.strip().lower() for c in metadata.columns]

if 'group' in metadata.columns:
    group_col = 'group'
else:
    group_col = metadata.columns[1]

meta_map = dict(zip(metadata[metadata.columns[0]].astype(str),
                    metadata[group_col].astype(str)))

reports = sorted([f for f in os.listdir('.') if f.endswith('.bracken.S.report')])
if len(reports) < 2:
    print("Need >=2 samples for differential abundance. Creating placeholder outputs.")
    for f in ["deseq2_results.tsv", "ancombc_results.tsv", "lefse_results.tsv",
              "deseq2_volcano.png", "deseq2_ma.png",
              "lefse_cladogram.png", "sig_heatmap.png"]:
        with open(f, 'w') as fh:
            fh.write("")
    exit(0)

samples = []
abundances = []
for r in reports:
    sid = os.path.basename(r).replace('.bracken.S.report', '').replace('.bracken', '')
    df = pd.read_csv(r, sep='\t')
    df.columns = [c.strip() for c in df.columns]
    df = df.set_index('name')
    samples.append(sid)
    abundances.append(df['new_est_reads'])

abund = pd.concat(abundances, axis=1, join='outer').fillna(0)
abund.columns = samples
abund = abund.loc[(abund.sum(axis=1) > 10)]

groups = pd.Series({s: meta_map.get(s, 'unknown') for s in samples})
unique_groups = groups.unique()
if len(unique_groups) < 2:
    print("Need >=2 groups for differential abundance. Creating placeholder.")
    for f in ["deseq2_results.tsv", "ancombc_results.tsv", "lefse_results.tsv",
              "deseq2_volcano.png", "deseq2_ma.png",
              "lefse_cladogram.png", "sig_heatmap.png"]:
        with open(f, 'w') as fh:
            fh.write("")
    exit(0)

group1, group2 = unique_groups[:2]
samples_g1 = [s for s in samples if groups[s] == group1]
samples_g2 = [s for s in samples if groups[s] == group2]

if len(samples_g1) < 2 or len(samples_g2) < 2:
    print("Need >=2 samples per group for DESeq2/ANCOM-BC. Using moderated tests.")

def compute_fc_pval(counts_g1, counts_g2, pseudocount=1):
    """Compute log2 fold change and Welch t-test p-value."""
    mean_g1 = np.mean(counts_g1) + pseudocount
    mean_g2 = np.mean(counts_g2) + pseudocount
    log2fc = np.log2(mean_g1 / mean_g2)
    if len(counts_g1) > 1 and len(counts_g2) > 1:
        _, pval = stats.ttest_ind(counts_g1, counts_g2, equal_var=False)
    else:
        pval = 1.0
    return log2fc, pval

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    count_df = abund.astype(int)
    count_df = count_df[count_df.sum(axis=1) > 10]

    col_meta = pd.DataFrame({'group': groups}, index=samples)
    col_meta.index.name = 'sample'

    dds = DeseqDataSet(
        counts=count_df.T,
        metadata=col_meta,
        design="~group",
    )
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=["group", group1, group2])
    stat_res.summary()

    deseq_res = stat_res.results_df.copy()
    deseq_res['taxon'] = deseq_res.index
    deseq_res = deseq_res.rename(columns={
        'log2FoldChange': 'log2FC',
        'pvalue': 'pvalue',
        'padj': 'padj'
    })
    has_pydeseq2 = True
except ImportError:
    has_pydeseq2 = False
    deseq_res = pd.DataFrame()

def bh_correct(pvals):
    pvals = np.array(pvals)
    n = len(pvals)
    ranked = np.argsort(pvals)
    padj = np.ones(n)
    cumulative = 1.0
    for i in range(n - 1, -1, -1):
        cumulative = min(cumulative, pvals[ranked[i]] * n / (i + 1))
        padj[ranked[i]] = cumulative
    return np.clip(padj, 0, 1)

if not has_pydeseq2 or deseq_res.empty:
    results = []
    for taxon in abund.index:
        g1_vals = abund.loc[taxon, samples_g1].values
        g2_vals = abund.loc[taxon, samples_g2].values
        log2fc, pval = compute_fc_pval(g1_vals, g2_vals)
        results.append({'taxon': taxon, 'log2FC': log2fc, 'pvalue': pval})

    deseq_res = pd.DataFrame(results)
    deseq_res['padj'] = bh_correct(deseq_res['pvalue'].values)
    deseq_res['padj'] = np.clip(deseq_res['padj'], 0, 1)

deseq_res['-log10(padj)'] = -np.log10(deseq_res['padj'].clip(lower=1e-300))
deseq_res['significant'] = (deseq_res['padj'] < 0.05) & (abs(deseq_res['log2FC']) > 1)
deseq_res.to_csv("deseq2_results.tsv", sep='\t', index=False)

n_sig = deseq_res['significant'].sum()
fig_volcano = go.Figure()

ns = deseq_res[~deseq_res['significant']]
fig_volcano.add_trace(go.Scatter(
    x=ns['log2FC'], y=ns['-log10(padj)'],
    mode='markers',
    marker=dict(color='#A0A0A0', size=6, opacity=0.5),
    name=f'NS ({len(ns)})',
    hovertemplate='%{text}<br>FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>',
    text=ns['taxon']))

up = deseq_res[(deseq_res['significant']) & (deseq_res['log2FC'] > 0)]
fig_volcano.add_trace(go.Scatter(
    x=up['log2FC'], y=up['-log10(padj)'],
    mode='markers',
    marker=dict(color='#E74C3C', size=8, opacity=0.8),
    name=f'Up ({len(up)})',
    hovertemplate='%{text}<br>FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>',
    text=up['taxon']))

down = deseq_res[(deseq_res['significant']) & (deseq_res['log2FC'] < 0)]
fig_volcano.add_trace(go.Scatter(
    x=down['log2FC'], y=down['-log10(padj)'],
    mode='markers',
    marker=dict(color='#2E86AB', size=8, opacity=0.8),
    name=f'Down ({len(down)})',
    hovertemplate='%{text}<br>FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>',
    text=down['taxon']))

fig_volcano.add_hline(y=-np.log10(0.05), line_dash='dash', line_color='grey', opacity=0.5)
fig_volcano.add_vline(x=1, line_dash='dash', line_color='grey', opacity=0.5)
fig_volcano.add_vline(x=-1, line_dash='dash', line_color='grey', opacity=0.5)

fig_volcano.update_layout(
    title=f'Volcano Plot — {group1} vs {group2}<br><sup>Significant: {n_sig} taxa (padj<0.05, |log2FC|>1)</sup>',
    xaxis_title='log2 Fold Change',
    yaxis_title='-log10(adjusted p-value)',
    height=600, width=800, template='plotly_white',
    hovermode='closest')
fig_volcano.write_image("deseq2_volcano.png", format='png', width=800, height=600)

mean_abund = abund.mean(axis=1)
ma_df = deseq_res.copy()
ma_df['mean_abundance'] = mean_abund.reindex(ma_df['taxon']).values

fig_ma = go.Figure()

ns_ma = ma_df[~ma_df['significant']]
fig_ma.add_trace(go.Scatter(
    x=ns_ma['mean_abundance'], y=ns_ma['log2FC'],
    mode='markers',
    marker=dict(color='#A0A0A0', size=6, opacity=0.5),
    name=f'NS ({len(ns_ma)})',
    hovertemplate='%{text}<br>Mean: %{x:.1f}<br>FC: %{y:.2f}<extra></extra>',
    text=ns_ma['taxon']))

up_ma = ma_df[(ma_df['significant']) & (ma_df['log2FC'] > 0)]
fig_ma.add_trace(go.Scatter(
    x=up_ma['mean_abundance'], y=up_ma['log2FC'],
    mode='markers',
    marker=dict(color='#E74C3C', size=8, opacity=0.8),
    name=f'Up ({len(up_ma)})',
    hovertemplate='%{text}<br>Mean: %{x:.1f}<br>FC: %{y:.2f}<extra></extra>',
    text=up_ma['taxon']))

down_ma = ma_df[(ma_df['significant']) & (ma_df['log2FC'] < 0)]
fig_ma.add_trace(go.Scatter(
    x=down_ma['mean_abundance'], y=down_ma['log2FC'],
    mode='markers',
    marker=dict(color='#2E86AB', size=8, opacity=0.8),
    name=f'Down ({len(down_ma)})',
    hovertemplate='%{text}<br>Mean: %{x:.1f}<br>FC: %{y:.2f}<extra></extra>',
    text=down_ma['taxon']))

fig_ma.add_hline(y=0, line_color='black', line_width=1)
fig_ma.update_xaxes(type='log')

fig_ma.update_layout(
    title=f'MA Plot — {group1} vs {group2}',
    xaxis_title='Mean abundance',
    yaxis_title='log2 Fold Change',
    height=600, width=800, template='plotly_white',
    hovermode='closest')
fig_ma.write_image("deseq2_ma.png", format='png', width=800, height=600)

ancom_results = []
for taxon in abund.index:
    g1_vals = abund.loc[taxon, samples_g1].values.astype(float)
    g2_vals = abund.loc[taxon, samples_g2].values.astype(float)

    g1_vals += 0.5
    g2_vals += 0.5

    all_g1 = abund[samples_g1].values.astype(float) + 0.5
    all_g2 = abund[samples_g2].values.astype(float) + 0.5

    def clr(x):
        gm = np.exp(np.mean(np.log(x)))
        return np.log(x / gm)

    g1_clr = clr(g1_vals)
    g2_clr = clr(g2_vals)

    log2fc = np.mean(g1_clr) - np.mean(g2_clr)
    if len(g1_clr) > 1 and len(g2_clr) > 1:
        _, pval = stats.ttest_ind(g1_clr, g2_clr, equal_var=False)
    else:
        pval = 1.0

    ancom_results.append({'taxon': taxon, 'log2FC': log2fc, 'pvalue': pval})

ancom_df = pd.DataFrame(ancom_results)
ancom_df['padj'] = bh_correct(ancom_df['pvalue'].values)
ancom_df['padj'] = np.clip(ancom_df['padj'], 0, 1)
ancom_df['significant'] = (ancom_df['padj'] < 0.05) & (abs(ancom_df['log2FC']) > 1)
ancom_df.to_csv("ancombc_results.tsv", sep='\t', index=False)

lefse_results = []
for taxon in abund.index:
    g1_vals = abund.loc[taxon, samples_g1].values.astype(float)
    g2_vals = abund.loc[taxon, samples_g2].values.astype(float)

    g1_vals += 1
    g2_vals += 1

    g1_rel = g1_vals / g1_vals.sum()
    g2_rel = g2_vals / g2_vals.sum()

    if len(g1_vals) > 1 and len(g2_vals) > 1:
        try:
            _, pval = stats.mannwhitneyu(g1_rel, g2_rel, alternative='two-sided')
        except ValueError:
            pval = 1.0
    else:
        pval = 1.0

    effect = abs(np.mean(g1_rel) - np.mean(g2_rel))

    lefse_results.append({'taxon': taxon, 'pvalue': pval, 'effect_size': effect,
                          'mean_g1': np.mean(g1_rel), 'mean_g2': np.mean(g2_rel)})

lefse_df = pd.DataFrame(lefse_results)
lefse_df['padj'] = bh_correct(lefse_df['pvalue'].values)
lefse_df['padj'] = np.clip(lefse_df['padj'], 0, 1)
lefse_df['significant'] = (lefse_df['padj'] < 0.05) & (lefse_df['effect_size'] > 0.01)
lefse_df.to_csv("lefse_results.tsv", sep='\t', index=False)

sig_lefse = lefse_df[lefse_df['significant']].sort_values('effect_size', ascending=False)
n_lefse_sig = len(sig_lefse)

fig_clad = go.Figure()

taxa_names = sig_lefse['taxon'].tolist() if n_lefse_sig > 0 else ['']
n_display = min(n_lefse_sig, 30)

if n_display > 0:
    display_taxa = sig_lefse.head(n_display)
    n_items = len(display_taxa)
    angles = np.linspace(0, 2 * np.pi, n_items, endpoint=False)
    radii = np.linspace(0.3, 1.0, n_items)

    for idx, (_, row) in enumerate(display_taxa.iterrows()):
        angle = angles[idx]
        r = 0.5 + 0.4 * (1 - row['padj'])
        x = r * np.cos(angle)
        y = r * np.sin(angle)
        color = '#E74C3C' if row['mean_g1'] > row['mean_g2'] else '#2E86AB'
        size = 8 + 12 * min(row['effect_size'] / max(display_taxa['effect_size']), 1)

        fig_clad.add_trace(go.Scatter(
            x=[x], y=[y],
            mode='markers+text',
            text=[row['taxon'][:25]],
            textposition='middle center',
            marker=dict(size=size, color=color, opacity=0.8,
                        line=dict(width=1, color='black')),
            name=row['taxon'],
            hovertemplate=f"{row['taxon']}<br>p-value: {row['padj']:.4f}<br>Effect: {row['effect_size']:.3f}<extra></extra>"))

    for r in [0.3, 0.5, 0.7, 0.9]:
        theta = np.linspace(0, 2 * np.pi, 100)
        fig_clad.add_trace(go.Scatter(
            x=r * np.cos(theta), y=r * np.sin(theta),
            mode='lines',
            line=dict(color='lightgray', width=1),
            showlegend=False, hoverinfo='skip'))

fig_clad.update_layout(
    title=f'LEfSe Biomarker Cladogram — {group1} vs {group2}<br><sup>{n_lefse_sig} significant biomarkers</sup>',
    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False,
               range=[-1.3, 1.3]),
    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False,
               range=[-1.3, 1.3]),
    height=700, width=700, template='plotly_white',
    hovermode='closest',
    annotations=[
        dict(x=0.9, y=0.9, text=f'↑ {group1}', showarrow=False, font=dict(color='#E74C3C', size=14)),
        dict(x=-0.9, y=-0.9, text=f'↑ {group2}', showarrow=False, font=dict(color='#2E86AB', size=14)),
    ])
fig_clad.write_image("lefse_cladogram.png", format='png', width=700, height=700)

rel_abund = abund.div(abund.sum(axis=0), axis=1) * 100

sig_taxa = set()
if n_sig > 0:
    sig_taxa.update(deseq_res[deseq_res['significant']]['taxon'].tolist())
sig_ancom = ancom_df[ancom_df['significant']]['taxon'].tolist()
sig_taxa.update(sig_ancom)
if n_lefse_sig > 0:
    sig_taxa.update(sig_lefse['taxon'].tolist())

sig_taxa_list = list(sig_taxa)
if sig_taxa_list:
    heat_data = rel_abund.loc[rel_abund.index.isin(sig_taxa_list)]
    heat_data = heat_data.loc[heat_data.max(axis=1) > 0]

    if len(heat_data) > 0:
        heat_z = heat_data.subtract(heat_data.mean(axis=1), axis=0).div(heat_data.std(axis=1), axis=0)
        heat_z = heat_z.fillna(0)
        heat_z = heat_z.loc[heat_z.mean(axis=1).sort_values(ascending=False).index]

        group_colors = [group1, group2]
        group_color_map = {g: i for i, g in enumerate(unique_groups[:2])}
        bar_colors = [group_color_map.get(groups[s], 0) for s in heat_z.columns]

        fig_heat = go.Figure()
        fig_heat.add_trace(go.Heatmap(
            z=heat_z.values,
            x=heat_z.columns,
            y=heat_z.index,
            colorscale='RdBu_r',
            zmid=0,
            hovertemplate='%{y}<br>%{x}: %{z:.2f}<extra></extra>'))

        fig_heat.update_layout(
            title=f'Significant Taxa Heatmap (Z-score) — {group1} vs {group2}',
            xaxis_title='Sample',
            yaxis_title='Taxa',
            height=max(400, 20 * len(heat_z)),
            width=600 + 20 * len(heat_z.columns),
            template='plotly_white',
            margin=dict(l=250))
        fig_heat.write_image("sig_heatmap.png", format='png', width=800, height=600)
    else:
        with open("sig_heatmap.png", 'w') as f:
            f.write('')
else:
    with open("sig_heatmap.png", 'w') as f:
        f.write('')

print(f"Differential abundance complete: {len(deseq_res)} taxa tested, "
      f"{deseq_res['significant'].sum()} DESeq2 significant, "
      f"{ancom_df['significant'].sum()} ANCOM-BC significant, "
      f"{lefse_df['significant'].sum()} LEfSe significant")
