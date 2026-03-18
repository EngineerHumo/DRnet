import csv
import math
import runpy

import pandas as pd
import gseapy as gp

from pipeline_utils import ensure_dirs, read_gmt, load_ensembl_symbol_mapping, normalize_ensembl_id, log_message

cfg = runpy.run_path('00_config.py')
RAW_DIR, PROC_DIR, RESULT_DIR = cfg['RAW_DIR'], cfg['PROC_DIR'], cfg['RESULT_DIR']
SEED = cfg.get('RANDOM_SEED', 202501)


def _run_preranked_gsea(rank_df, gene_sets):
    rnk = pd.DataFrame({'gene': rank_df['gene'], 'score': rank_df['score']})
    pre = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        processes=1,
        permutation_num=2000,
        min_size=5,
        max_size=2000,
        seed=SEED,
        outdir=None,
        verbose=False,
    )
    tab = pre.res2d.reset_index().rename(columns={'Term': 'pathway'})

    rows = []
    for _, r in tab.iterrows():
        rows.append({
            'pathway': str(r['pathway']),
            'ES': float(r.get('ES', 0.0)),
            'NES': float(r.get('NES', 0.0)),
            'pvalue': float(r.get('NOM p-val', 1.0)),
            'padj': float(r.get('FDR q-val', 1.0)),
            'hit_genes': len(set(gene_sets.get(str(r['pathway']), [])).intersection(rank_df['gene'])),
            'geneset_size': len(set(gene_sets.get(str(r['pathway']), []))),
            'engine': 'gseapy_prerank',
        })
    return rows


def main():
    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'logs')

    with open(RESULT_DIR / 'tables' / 'lasso_selected_genes.csv', encoding='utf-8') as f:
        genes = [r.get('gene_symbol', r.get('gene', '')) for r in csv.DictReader(f)]
        genes = [g for g in genes if g]

    with open(RESULT_DIR / 'tables' / 'deg_primary_healthy_vs_npdr_pdr_dme.csv', encoding='utf-8') as f:
        deg = list(csv.DictReader(f))
    mapping = load_ensembl_symbol_mapping(RAW_DIR, PROC_DIR, ensembl_ids=[r['gene'] for r in deg])

    gene_score = {}
    for r in deg:
        eid = normalize_ensembl_id(r['gene'])
        sym = mapping.get(eid)
        if not sym:
            continue
        p = max(float(r['pvalue']), 1e-300)
        score = float(r['log2FC']) * (-math.log10(p))
        if sym not in gene_score or abs(score) > abs(gene_score[sym]):
            gene_score[sym] = score

    ranked = sorted(gene_score.items(), key=lambda x: x[1], reverse=True)[:6000]
    rank_df = {'gene': [g for g, _ in ranked], 'score': [float(s) for _, s in ranked]}

    hall = read_gmt(RAW_DIR / 'hallmark_human_gene_symbols.gmt')
    rows = _run_preranked_gsea(rank_df, hall)
    rows.sort(key=lambda x: (x['padj'], -abs(x['NES'])))

    with open(RESULT_DIR / 'tables' / 'gsea_summary_matrix.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['pathway', 'ES', 'NES', 'pvalue', 'padj', 'hit_genes', 'geneset_size', 'engine'])
        w.writeheader()
        w.writerows(rows)

    with open(RESULT_DIR / 'tables' / 'pathway_recurrence_summary.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['pathway', 'hit_genes', 'padj', 'NES'])
        w.writeheader()
        w.writerows([{'pathway': r['pathway'], 'hit_genes': r['hit_genes'], 'padj': r['padj'], 'NES': r['NES']} for r in rows])

    hall_sets = {k: set(v) for k, v in hall.items()}
    for g in genes:
        out = []
        for r in rows:
            if g in hall_sets.get(r['pathway'], set()):
                out.append({'gene_symbol': g, 'pathway': r['pathway'], 'NES': r['NES'], 'padj': r['padj']})
        with open(RESULT_DIR / 'tables' / f'gsea_{g}.csv', 'w', newline='', encoding='utf-8') as f:
            fields = ['gene_symbol', 'pathway', 'NES', 'padj']
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(out)

    with open(RESULT_DIR / 'tables' / 'ipa_input_selected_genes.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['gene_symbol'])
        w.writeheader()
        for g in genes:
            w.writerow({'gene_symbol': g})

    log_message('08_enrichment_analysis', f'ranked_genes={len(ranked)} pathways={len(rows)} selected_genes={len(genes)} engine=gseapy_prerank')


if __name__ == '__main__':
    main()
