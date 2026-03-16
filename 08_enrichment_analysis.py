import csv, runpy
from pipeline_utils import ensure_dirs, read_gmt, log_message

cfg = runpy.run_path('00_config.py')
RAW_DIR, RESULT_DIR = cfg['RAW_DIR'], cfg['RESULT_DIR']


def main():
    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'logs')
    with open(RESULT_DIR / 'tables' / 'lasso_selected_genes.csv', encoding='utf-8') as f:
        genes = [r.get('gene_symbol', r.get('gene', '')) for r in csv.DictReader(f)]
        genes = [g for g in genes if g]

    hall = read_gmt(RAW_DIR / 'hallmark_human_gene_symbols.gmt')
    summ = []
    for g in genes:
        rows = []
        for p, gl in hall.items():
            overlap = int(g in set(gl))
            rows.append({'gene_symbol': g, 'pathway': p, 'overlap': overlap, 'spearman_rank_score': overlap})
        with open(RESULT_DIR / 'tables' / f'gsea_{g}.csv', 'w', newline='', encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=rows[0].keys())
            w.writeheader()
            w.writerows(rows)
        summ.extend(rows)

    with open(RESULT_DIR / 'tables' / 'gsea_summary_matrix.csv', 'w', newline='', encoding='utf-8') as f:
        fields = summ[0].keys() if summ else ['gene_symbol', 'pathway', 'overlap', 'spearman_rank_score']
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(summ)

    rec = {}
    for r in summ:
        rec[r['pathway']] = rec.get(r['pathway'], 0) + r['overlap']
    with open(RESULT_DIR / 'tables' / 'pathway_recurrence_summary.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['pathway', 'hit_genes'])
        w.writeheader()
        w.writerows([{'pathway': k, 'hit_genes': v} for k, v in sorted(rec.items(), key=lambda x: -x[1])])

    with open(RESULT_DIR / 'tables' / 'ipa_input_selected_genes.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['gene_symbol'])
        w.writeheader()
        w.writerows([{'gene_symbol': g} for g in genes])

    log_message('08_enrichment_analysis', f'genes={len(genes)}')


if __name__ == '__main__':
    main()
