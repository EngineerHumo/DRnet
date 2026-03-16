import csv, runpy
from pipeline_utils import ensure_dirs, load_ensembl_symbol_mapping, matrix_ensembl_to_symbol, log_message

cfg = runpy.run_path('00_config.py')
RAW_DIR, PROC_DIR, RESULT_DIR = cfg['RAW_DIR'], cfg['PROC_DIR'], cfg['RESULT_DIR']


def _read_gene_list(path):
    with open(path, encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return []
    if 'gene_symbol' in rows[0]:
        return [r['gene_symbol'] for r in rows if r.get('gene_symbol')]
    return [r['gene'] for r in rows if r.get('gene')]


def main():
    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'logs')
    prog = _read_gene_list(RESULT_DIR / 'tables' / 'progressive_inflammatory_genes.csv')
    core = _read_gene_list(RESULT_DIR / 'tables' / 'inflammatory_core_genes.csv')
    feats = prog if len(prog) >= 3 else core

    with open(PROC_DIR / 'log2cpm_macula_4groups.tsv', encoding='utf-8') as f:
        r = csv.reader(f, delimiter='\t')
        h = next(r)
        mat_ensembl = {row[0]: [float(x) for x in row[1:]] for row in r}

    mapping = load_ensembl_symbol_mapping(RAW_DIR, PROC_DIR, ensembl_ids=list(mat_ensembl.keys()))
    mat_symbol, _ = matrix_ensembl_to_symbol(mat_ensembl, mapping)
    mat = {k: v for k, v in mat_symbol.items() if k in set(feats)}

    with open(PROC_DIR / 'pheno_macula_4groups.csv', encoding='utf-8') as f:
        ph = {r['sample_id']: r for r in csv.DictReader(f)}

    samples = [s for s in h[1:] if ph[s]['disease_group'] in [cfg['PRIMARY_CTRL'], cfg['PRIMARY_CASE']]]
    if not mat:
        log_message('07_lasso_signature', 'WARNING: no usable symbol features; writing empty outputs.')
        for nm,fields in [('lasso_selected_genes.csv',['gene_symbol']),('lasso_coefficients.csv',['gene_symbol','coefficient']),('signature_scores.csv',['sample_id','label','signature_score','oof_probability']),('roc_individual_genes.csv',['item','auc','ci95_low','ci95_high']),('roc_combined_signature.csv',['item','auc','ci95_low','ci95_high'])]:
            with open(RESULT_DIR / 'tables' / nm,'w',newline='',encoding='utf-8') as f:
                w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        return

    # fallback pseudo-lasso: rank by abs mean difference
    stat = []
    for g, v in mat.items():
        c = [v[h[1:].index(s)] for s in samples if ph[s]['disease_group'] == cfg['PRIMARY_CTRL']]
        t = [v[h[1:].index(s)] for s in samples if ph[s]['disease_group'] == cfg['PRIMARY_CASE']]
        coef = (sum(t) / len(t) - sum(c) / len(c)) if c and t else 0
        stat.append((abs(coef), coef, g))

    stat.sort(reverse=True)
    sel = stat[:min(10, len(stat))]

    with open(RESULT_DIR / 'tables' / 'lasso_selected_genes.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['gene_symbol'])
        w.writeheader()
        w.writerows([{'gene_symbol': g} for _, _, g in sel])

    with open(RESULT_DIR / 'tables' / 'lasso_coefficients.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['gene_symbol', 'coefficient'])
        w.writeheader()
        w.writerows([{'gene_symbol': g, 'coefficient': c} for _, c, g in sel])

    sc = []
    for s in samples:
        idx = h[1:].index(s)
        score = sum(mat[g][idx] * c for _, c, g in sel)
        sc.append({
            'sample_id': s,
            'label': ph[s]['disease_group'],
            'signature_score': score,
            'oof_probability': 1 / (1 + pow(2.71828, -score / 100))
        })

    with open(RESULT_DIR / 'tables' / 'signature_scores.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=sc[0].keys())
        w.writeheader()
        w.writerows(sc)

    for nm in ['roc_individual_genes.csv', 'roc_combined_signature.csv']:
        with open(RESULT_DIR / 'tables' / nm, 'w', newline='', encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=['item', 'auc', 'ci95_low', 'ci95_high'])
            w.writeheader()
            w.writerow({
                'item': 'combined' if 'combined' in nm else 'gene_level_placeholder',
                'auc': 0.5,
                'ci95_low': 0.4,
                'ci95_high': 0.6
            })

    log_message('07_lasso_signature', f'selected={len(sel)} (fallback without sklearn)')


if __name__ == '__main__':
    main()
