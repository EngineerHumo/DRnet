import csv, runpy
from pipeline_utils import (
    ensure_dirs, read_gmt, ssgsea_score, mann_whitney_u, spearman,
    load_ensembl_symbol_mapping, matrix_ensembl_to_symbol, log_message
)

cfg = runpy.run_path('00_config.py')
RAW_DIR, PROC_DIR, RESULT_DIR = cfg['RAW_DIR'], cfg['PROC_DIR'], cfg['RESULT_DIR']


def main():
    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'logs')
    hallmark_genes = set(read_gmt(RAW_DIR / 'hallmark_inflammatory_response.gmt')['HALLMARK_INFLAMMATORY_RESPONSE'])

    with open(PROC_DIR / 'log2cpm_macula_4groups.tsv', encoding='utf-8') as f:
        r = csv.reader(f, delimiter='\t')
        h = next(r)
        mat_ensembl = {row[0]: [float(x) for x in row[1:]] for row in r}

    mapping = load_ensembl_symbol_mapping(RAW_DIR, PROC_DIR, ensembl_ids=list(mat_ensembl.keys()))
    mat_symbol, dedup_log = matrix_ensembl_to_symbol(mat_ensembl, mapping)

    overlap = hallmark_genes.intersection(set(mat_symbol.keys()))
    with open(RESULT_DIR / 'logs' / '05_inflammation_scoring_mapping_check.md', 'w', encoding='utf-8') as f:
        f.write(
            f"# Mapping check\n\n"
            f"- ensembl rows: {len(mat_ensembl)}\n"
            f"- mapped symbol rows: {len(mat_symbol)}\n"
            f"- hallmark genes: {len(hallmark_genes)}\n"
            f"- overlap genes used in ssGSEA: {len(overlap)}\n"
            f"- deduplicated symbols: {len(dedup_log)}\n"
        )

    if len(overlap) == 0:
        log_message('05_inflammation_scoring', 'WARNING: no symbol overlap; scores will be set to 0. Provide data_raw/ensembl_to_symbol.csv for valid biology.')

    with open(PROC_DIR / 'pheno_macula_4groups.csv', encoding='utf-8') as f:
        ph = {r['sample_id']: r for r in csv.DictReader(f)}

    rows = []
    for i, s in enumerate(h[1:]):
        expr = {g: v[i] for g, v in mat_symbol.items()}
        rows.append({
            'sample_id': s,
            'group': ph[s]['disease_group'],
            'severity_code': ph[s]['severity_code'],
            'inflammation_ssgsea': ssgsea_score(expr, hallmark_genes)
        })

    with open(RESULT_DIR / 'tables' / 'inflammation_ssgsea_scores.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    x = [r['inflammation_ssgsea'] for r in rows if r['group'] == cfg['PRIMARY_CTRL']]
    y = [r['inflammation_ssgsea'] for r in rows if r['group'] == cfg['PRIMARY_CASE']]
    _, p = mann_whitney_u(x, y)
    with open(RESULT_DIR / 'tables' / 'inflammation_group_comparison.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['contrast', 'pvalue'])
        w.writeheader()
        w.writerow({'contrast': 'healthy control vs NPDR/PDR + DME', 'pvalue': p})

    rho, p2 = spearman([int(r['severity_code']) for r in rows], [r['inflammation_ssgsea'] for r in rows])
    with open(RESULT_DIR / 'tables' / 'inflammation_severity_trend.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['rho', 'pvalue'])
        w.writeheader()
        w.writerow({'rho': rho, 'pvalue': p2})

    log_message('05_inflammation_scoring', f'done, overlap_genes={len(overlap)}')


if __name__ == '__main__':
    main()
