import csv
import runpy

from pipeline_utils import ensure_dirs, load_ensembl_symbol_mapping, matrix_ensembl_to_symbol, log_message

cfg = runpy.run_path('00_config.py')
RAW_DIR, PROC_DIR, RESULT_DIR = cfg['RAW_DIR'], cfg['PROC_DIR'], cfg['RESULT_DIR']
SEED = cfg.get('RANDOM_SEED', 202501)


def _read_gene_list(path):
    with open(path, encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return []
    if 'gene_symbol' in rows[0]:
        return [r['gene_symbol'] for r in rows if r.get('gene_symbol')]
    return [r['gene'] for r in rows if r.get('gene')]


def _bootstrap_auc_ci(y_true, y_prob, n_boot=500):
    import numpy as np
    from sklearn.metrics import roc_auc_score

    rng = np.random.default_rng(SEED)
    y_true = np.asarray(y_true, dtype=int)
    y_prob = np.asarray(y_prob, dtype=float)
    idx = np.arange(len(y_true))
    vals = []
    for _ in range(n_boot):
        bs = rng.choice(idx, size=len(idx), replace=True)
        ys = y_true[bs]
        if len(np.unique(ys)) < 2:
            continue
        vals.append(roc_auc_score(ys, y_prob[bs]))
    if not vals:
        return 0.5, 0.5
    vals.sort()
    return float(vals[int(0.025 * len(vals))]), float(vals[int(0.975 * len(vals)) - 1])


def main():
    try:
        import numpy as np
        from sklearn.linear_model import LogisticRegressionCV
        from sklearn.metrics import roc_auc_score
        from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
    except ImportError as e:
        raise ImportError('07_lasso_signature.py 需要 numpy + scikit-learn。') from e

    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'logs')
    prog = _read_gene_list(RESULT_DIR / 'tables' / 'progressive_inflammatory_genes.csv')
    core = _read_gene_list(RESULT_DIR / 'tables' / 'inflammatory_core_genes.csv')
    feats = sorted(set(prog if len(prog) >= 6 else core))

    with open(PROC_DIR / 'log2cpm_macula_4groups.tsv', encoding='utf-8') as f:
        r = csv.reader(f, delimiter='\t')
        h = next(r)
        mat_ensembl = {row[0]: [float(x) for x in row[1:]] for row in r}

    mapping = load_ensembl_symbol_mapping(RAW_DIR, PROC_DIR, ensembl_ids=list(mat_ensembl.keys()))
    mat_symbol, _ = matrix_ensembl_to_symbol(mat_ensembl, mapping)

    with open(PROC_DIR / 'pheno_macula_4groups.csv', encoding='utf-8') as f:
        ph = {r['sample_id']: r for r in csv.DictReader(f)}

    sample_ids = [s for s in h[1:] if ph[s]['disease_group'] in [cfg['PRIMARY_CTRL'], cfg['PRIMARY_CASE']]]
    y = np.asarray([1 if ph[s]['disease_group'] == cfg['PRIMARY_CASE'] else 0 for s in sample_ids], dtype=int)

    available = [g for g in feats if g in mat_symbol]
    if len(available) < 3:
        with open(RESULT_DIR / 'tables' / 'deg_primary_healthy_vs_npdr_pdr_dme.csv', encoding='utf-8') as f:
            top_ens = [r['gene'] for r in list(csv.DictReader(f))[:80]]
        ordered, seen = [], set()
        for eid in top_ens:
            sym = mapping.get(eid.split('.')[0])
            if sym and sym in mat_symbol and sym not in seen:
                ordered.append(sym)
                seen.add(sym)
            if len(ordered) >= 30:
                break
        available = ordered

    sample_idx = {s: i for i, s in enumerate(h[1:])}
    X = np.asarray([[mat_symbol[g][sample_idx[s]] for g in available] for s in sample_ids], dtype=float)

    cs = np.logspace(-3, 2, 30)
    inner_cv = RepeatedStratifiedKFold(n_splits=4, n_repeats=5, random_state=SEED)
    outer_cv = StratifiedKFold(n_splits=min(5, max(2, len(y) // 4)), shuffle=True, random_state=SEED)

    oof = np.full(len(y), 0.5, dtype=float)
    fold_c = []
    for tr_idx, te_idx in outer_cv.split(X, y):
        model = Pipeline([
            ('scaler', StandardScaler()),
            ('clf', LogisticRegressionCV(
                Cs=cs,
                cv=inner_cv,
                penalty='l1',
                solver='saga',
                scoring='roc_auc',
                max_iter=3000,
                n_jobs=None,
                random_state=SEED,
                refit=True,
            ))
        ])
        model.fit(X[tr_idx], y[tr_idx])
        prob = model.predict_proba(X[te_idx])[:, 1]
        oof[te_idx] = prob
        fold_c.append(float(model.named_steps['clf'].C_[0]))

    best_c = float(np.median(fold_c)) if fold_c else float(cs[len(cs) // 2])
    final_model = Pipeline([
        ('scaler', StandardScaler()),
        ('clf', LogisticRegressionCV(
            Cs=[best_c],
            cv=inner_cv,
            penalty='l1',
            solver='saga',
            scoring='roc_auc',
            max_iter=3000,
            n_jobs=None,
            random_state=SEED,
            refit=True,
        ))
    ])
    final_model.fit(X, y)

    coefs = final_model.named_steps['clf'].coef_[0]
    selected = [(g, float(c)) for g, c in zip(available, coefs) if abs(c) > 1e-10]
    if not selected:
        selected = sorted([(g, float(c)) for g, c in zip(available, coefs)], key=lambda x: -abs(x[1]))[:10]

    with open(RESULT_DIR / 'tables' / 'lasso_selected_genes.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['gene_symbol'])
        w.writeheader()
        for g, _ in sorted(selected, key=lambda x: -abs(x[1])):
            w.writerow({'gene_symbol': g})

    with open(RESULT_DIR / 'tables' / 'lasso_coefficients.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['gene_symbol', 'coefficient'])
        w.writeheader()
        for g, c in sorted(selected, key=lambda x: -abs(x[1])):
            w.writerow({'gene_symbol': g, 'coefficient': c})

    full_prob = final_model.predict_proba(X)[:, 1]
    decision = final_model.decision_function(X)
    scores = [
        {
            'sample_id': s,
            'label': ph[s]['disease_group'],
            'signature_score': float(decision[i]),
            'oof_probability': float(oof[i]),
            'full_probability': float(full_prob[i]),
        }
        for i, s in enumerate(sample_ids)
    ]
    with open(RESULT_DIR / 'tables' / 'signature_scores.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=scores[0].keys())
        w.writeheader()
        w.writerows(scores)

    scaler = final_model.named_steps['scaler']
    Xz = scaler.transform(X)
    roc_rows = []
    for g, _ in selected:
        j = available.index(g)
        gene_prob = 1.0 / (1.0 + np.exp(-Xz[:, j]))
        auc = float(roc_auc_score(y, gene_prob))
        lo, hi = _bootstrap_auc_ci(y, gene_prob)
        roc_rows.append({'item': g, 'auc': auc, 'ci95_low': lo, 'ci95_high': hi})

    with open(RESULT_DIR / 'tables' / 'roc_individual_genes.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['item', 'auc', 'ci95_low', 'ci95_high'])
        w.writeheader()
        w.writerows(sorted(roc_rows, key=lambda x: -x['auc']))

    auc = float(roc_auc_score(y, oof))
    lo, hi = _bootstrap_auc_ci(y, oof)
    with open(RESULT_DIR / 'tables' / 'roc_combined_signature.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['item', 'auc', 'ci95_low', 'ci95_high', 'lambda'])
        w.writeheader()
        w.writerow({'item': 'combined_signature_oof', 'auc': auc, 'ci95_low': lo, 'ci95_high': hi, 'lambda': 1.0 / best_c})

    log_message('07_lasso_signature', f'features={len(available)} selected={len(selected)} auc={auc:.3f} C={best_c:.6g}')


if __name__ == '__main__':
    main()
