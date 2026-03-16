import csv
import runpy
from collections import Counter, defaultdict
from pipeline_utils import ensure_dirs, log_message
from simple_plot import Canvas, scale

cfg = runpy.run_path('00_config.py')
RESULT_DIR = cfg['RESULT_DIR']
PROC_DIR = cfg['PROC_DIR']


def read_csv(path):
    with open(path, encoding='utf-8') as f:
        return list(csv.DictReader(f))


def save(canvas, name):
    canvas.save_png(RESULT_DIR / 'figures' / f'{name}.png')
    canvas.save_pdf(RESULT_DIR / 'figures' / f'{name}.pdf')


def draw_bar_from_pairs(pairs, name, color=(80, 145, 220)):
    c = Canvas(1100, 700)
    c.line(90, 620, 1030, 620, (0, 0, 0), 2)
    c.line(90, 620, 90, 90, (0, 0, 0), 2)
    vals = [v for _, v in pairs] or [1]
    vmax = max(vals)
    bw = max(8, int(880 / max(1, len(pairs))))
    for i, (_, v) in enumerate(pairs):
        x = 110 + i * bw
        y = int(scale(v, 0, vmax, 620, 120))
        c.fill_rect(x, y, x + bw - 2, 620, color)
    save(c, name)


def draw_scatter(rows, xk, yk, groupk, name, palette=None):
    c = Canvas(1100, 700)
    default_palette = {'healthy control': (65, 105, 225), 'diabetic': (34, 139, 34), 'NPDR': (255, 140, 0), 'NPDR/PDR + DME': (220, 20, 60)}
    palette = palette or default_palette
    xs = [float(r[xk]) for r in rows]
    ys = [float(r[yk]) for r in rows]
    c.line(90, 620, 1030, 620, (0, 0, 0), 2)
    c.line(90, 620, 90, 90, (0, 0, 0), 2)
    for r in rows:
        x = int(scale(float(r[xk]), min(xs), max(xs), 120, 1000))
        y = int(scale(float(r[yk]), min(ys), max(ys), 600, 120))
        grp = str(r.get(groupk, ''))
        if grp in palette:
            color = palette[grp]
        else:
            # deterministic fallback color for unseen classes (avoid all-gray collapse)
            h = sum(ord(ch) for ch in grp)
            color = (40 + (h * 53) % 180, 40 + (h * 97) % 180, 40 + (h * 29) % 180)
        c.circle(x, y, 5, color)
    save(c, name)


def draw_heatmap(rows, row_key, col_key, val_key, name):
    rnames = sorted({r[row_key] for r in rows})
    cnames = sorted({r[col_key] for r in rows})
    tab = {(r[row_key], r[col_key]): float(r[val_key]) for r in rows}
    vals = [v for v in tab.values()] or [0.0]
    vmin, vmax = min(vals), max(vals)
    c = Canvas(1000, 800)
    x0, y0, w, h = 80, 80, 840, 640
    cw = max(1, w // max(1, len(cnames)))
    ch = max(1, h // max(1, len(rnames)))
    for i, rn in enumerate(rnames):
        for j, cn in enumerate(cnames):
            v = tab.get((rn, cn), 0.0)
            t = 0.5 if vmax == vmin else (v - vmin) / (vmax - vmin)
            col = (int(255 * (1 - t)), int(255 * (1 - t)), 255)
            c.fill_rect(x0 + j * cw, y0 + i * ch, x0 + (j + 1) * cw - 1, y0 + (i + 1) * ch - 1, col)
    save(c, name)


def main():
    ensure_dirs(RESULT_DIR / 'figures', RESULT_DIR / 'logs')

    # Figure1
    with open(PROC_DIR / 'pheno_macula_4groups.csv', encoding='utf-8') as f:
        ph = list(csv.DictReader(f))
    cnt = Counter(r['disease_group'] for r in ph)
    draw_bar_from_pairs(sorted(cnt.items()), 'Figure1_sample_composition', (90, 150, 230))
    pca = read_csv(RESULT_DIR / 'tables' / 'pca_coordinates.csv')
    draw_scatter(pca, 'PC1', 'PC2', 'group', 'Figure1_pca_macula_4groups')

    ms = read_csv(RESULT_DIR / 'tables' / 'master_analysis_summary.csv') if (RESULT_DIR / 'tables' / 'master_analysis_summary.csv').exists() else []
    if ms:
        draw_bar_from_pairs([(r['table'], float(r['rows'])) for r in ms[:30]], 'Figure1_workflow', (160, 160, 160))
    else:
        draw_bar_from_pairs([('steps', 12)], 'Figure1_workflow', (160, 160, 160))

    # Figure2
    infl = read_csv(RESULT_DIR / 'tables' / 'inflammation_ssgsea_scores.csv')
    grp = defaultdict(list)
    for r in infl:
        grp[r['group']].append(float(r['inflammation_ssgsea']))
    draw_bar_from_pairs([(k, sum(v) / len(v)) for k, v in sorted(grp.items())], 'Figure2_inflammation_ssgsea_4groups', (220, 140, 90))

    deg_all = read_csv(RESULT_DIR / 'tables' / 'deg_primary_healthy_vs_npdr_pdr_dme.csv')
    deg_plot = deg_all[:1000]
    volc = [{'x': float(r['log2FC']), 'y': -__import__('math').log10(max(float(r['pvalue']), 1e-300)), 'group': 'sig' if int(r['significant']) else 'ns'} for r in deg_plot]
    draw_scatter(volc, 'x', 'y', 'group', 'Figure2_primary_volcano', palette={'sig': (220, 20, 60), 'ns': (160, 160, 160)})

    top_deg = deg_all[:40]
    heat_rows = [{'gene': r['gene'], 'metric': 'log2FC', 'val': float(r['log2FC'])} for r in top_deg]
    draw_heatmap(heat_rows, 'gene', 'metric', 'val', 'Figure2_deg_inflammatory_heatmap')

    # overlap counts
    with open(RESULT_DIR / 'tables' / 'inflammatory_core_genes.csv', encoding='utf-8') as f:
        core = [r for r in csv.DictReader(f)]
    draw_bar_from_pairs([('primary_deg_sig', sum(int(r['significant']) for r in deg_all)), ('inflam_core', len(core))], 'Figure2_overlap_primarydeg_inflammation', (150, 200, 120))

    # Figure3
    coef = read_csv(RESULT_DIR / 'tables' / 'lasso_coefficients.csv')
    draw_bar_from_pairs([(r['gene_symbol'], abs(float(r['coefficient']))) for r in coef[:25]], 'Figure3_lasso_path', (200, 100, 120))
    roc = read_csv(RESULT_DIR / 'tables' / 'roc_combined_signature.csv')
    draw_bar_from_pairs([(r['item'], float(r['auc'])) for r in roc], 'Figure3_cv_curve', (120, 120, 220))

    sig_scores = read_csv(RESULT_DIR / 'tables' / 'signature_scores.csv')
    draw_scatter(sig_scores, 'signature_score', 'oof_probability', 'label', 'Figure3_selected_gene_expression_4groups')

    ind_roc = read_csv(RESULT_DIR / 'tables' / 'roc_individual_genes.csv')
    draw_bar_from_pairs([(r['item'], float(r['auc'])) for r in ind_roc[:20]], 'Figure3_individual_roc', (120, 180, 140))
    draw_bar_from_pairs([(r['item'], float(r['auc'])) for r in roc], 'Figure3_combined_roc', (80, 120, 200))

    # Figure4
    gsea = read_csv(RESULT_DIR / 'tables' / 'gsea_summary_matrix.csv')
    draw_bar_from_pairs([(r['pathway'], abs(float(r['NES']))) for r in gsea[:20]], 'Figure4_recurrent_pathway_dotplot', (180, 100, 200))
    draw_heatmap([{'pathway': r['pathway'], 'metric': 'NES', 'val': float(r['NES'])} for r in gsea[:40]], 'pathway', 'metric', 'val', 'Figure4_selected_gene_gsea_heatmap')
    draw_bar_from_pairs([(r['pathway'], float(r['hit_genes'])) for r in gsea[:20]], 'Figure4_gene_pathway_network_optional', (200, 160, 80))

    # Figure5
    imm = read_csv(RESULT_DIR / 'tables' / 'immune_primary_comparison.csv')
    draw_bar_from_pairs([(r['cell_type'], -__import__('math').log10(max(float(r['padj']), 1e-300))) for r in imm[:28]], 'Figure5_immune_score_comparison', (100, 190, 190))

    cor = read_csv(RESULT_DIR / 'tables' / 'gene_immune_correlations.csv')
    if cor:
        draw_heatmap([{'gene': r['gene_symbol'], 'cell': r['cell_type'], 'val': float(r['rho'])} for r in cor], 'gene', 'cell', 'val', 'Figure5_gene_immune_corr_heatmap')
    else:
        draw_bar_from_pairs([('no_sig_celltype', 0)], 'Figure5_gene_immune_corr_heatmap')

    # Supplementary
    qc = read_csv(RESULT_DIR / 'tables' / 'qc_metrics.csv')
    draw_bar_from_pairs([(r['sample_id'], float(r['library_size'])) for r in qc], 'Supp_QC', (120, 150, 220))
    draw_heatmap([{'sample': r['sample_id'], 'metric': 'expr_median', 'val': float(r['expr_median'])} for r in qc], 'sample', 'metric', 'val', 'Supp_clustering')
    draw_heatmap([{'gene': r['gene'], 'metric': 'padj', 'val': -__import__('math').log10(max(float(r['padj']), 1e-300))} for r in deg_all[:50]], 'gene', 'metric', 'val', 'Supp_expanded_deg_heatmap')
    draw_bar_from_pairs([(r['gene_symbol'], 1) for r in core[:30]], 'Supp_all_inflammatory_core_heatmap', (200, 80, 80))
    draw_bar_from_pairs([(r['pathway'], abs(float(r['NES']))) for r in gsea[:30]], 'Supp_single_gene_gsea', (140, 120, 200))
    if cor:
        draw_scatter([{'x': float(r['rho']), 'y': -__import__('math').log10(max(float(r['pvalue']), 1e-300)), 'grp': r['gene_symbol']} for r in cor], 'x', 'y', 'grp', 'Supp_gene_immune_scatter')
    else:
        draw_bar_from_pairs([('empty', 0)], 'Supp_gene_immune_scatter')

    log_message('10_make_figures', 'generated data-driven figures from tables')


if __name__ == '__main__':
    main()
