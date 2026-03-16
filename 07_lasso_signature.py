import csv
import math
import random
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


def auc_score(y, p):
    pos = [pp for pp, yy in zip(p, y) if yy == 1]
    neg = [pp for pp, yy in zip(p, y) if yy == 0]
    if not pos or not neg:
        return 0.5
    wins = 0.0
    for pp in pos:
        for pn in neg:
            if pp > pn:
                wins += 1
            elif pp == pn:
                wins += 0.5
    return wins / (len(pos) * len(neg))


def bootstrap_auc_ci(y, p, n_boot=300):
    rnd = random.Random(SEED)
    idx = list(range(len(y)))
    vals = []
    for _ in range(n_boot):
        bs = [rnd.choice(idx) for _ in idx]
        vals.append(auc_score([y[i] for i in bs], [p[i] for i in bs]))
    vals.sort()
    return vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals)) - 1]


def fit_scaler(X):
    n = len(X)
    p = len(X[0])
    means, stds = [0.0] * p, [1.0] * p
    for j in range(p):
        col = [X[i][j] for i in range(n)]
        m = sum(col) / n
        v = sum((a - m) ** 2 for a in col) / max(1, n - 1)
        means[j] = m
        stds[j] = math.sqrt(v) if v > 1e-12 else 1.0
    return means, stds


def transform_scaler(X, means, stds):
    return [[(row[j] - means[j]) / stds[j] for j in range(len(row))] for row in X]


def fit_l1_logistic(X, y, lam=0.1, lr=0.03, max_iter=800):
    n, p = len(X), len(X[0])
    w = [0.0] * p
    b = 0.0
    for _ in range(max_iter):
        grad_w = [0.0] * p
        grad_b = 0.0
        for i in range(n):
            z = b + sum(w[j] * X[i][j] for j in range(p))
            pr = 1.0 / (1.0 + math.exp(-max(min(z, 30), -30)))
            d = pr - y[i]
            grad_b += d
            for j in range(p):
                grad_w[j] += d * X[i][j]
        grad_b /= n
        for j in range(p):
            grad_w[j] /= n
        b -= lr * grad_b
        for j in range(p):
            v = w[j] - lr * grad_w[j]
            if v > lr * lam:
                w[j] = v - lr * lam
            elif v < -lr * lam:
                w[j] = v + lr * lam
            else:
                w[j] = 0.0
    return w, b


def predict_prob(X, w, b):
    out = []
    for row in X:
        z = b + sum(w[j] * row[j] for j in range(len(w)))
        z = max(min(z, 30), -30)
        out.append(1.0 / (1.0 + math.exp(-z)))
    return out


def kfold_indices(n, k=5):
    idx = list(range(n))
    rnd = random.Random(SEED)
    rnd.shuffle(idx)
    folds = [[] for _ in range(k)]
    for i, v in enumerate(idx):
        folds[i % k].append(v)
    return folds


def choose_lambda_inner_cv(X_train_raw, y_train, lambda_grid):
    n = len(X_train_raw)
    if n < 6:
        # small-n fallback: fixed split for held-out tuning
        split = max(1, n // 3)
        tr_idx = list(range(split, n))
        va_idx = list(range(0, split))
        if len(tr_idx) < 2 or len(va_idx) < 1:
            return lambda_grid[len(lambda_grid) // 2]

        Xtr_raw = [X_train_raw[i] for i in tr_idx]
        ytr = [y_train[i] for i in tr_idx]
        Xva_raw = [X_train_raw[i] for i in va_idx]
        yva = [y_train[i] for i in va_idx]

        means, stds = fit_scaler(Xtr_raw)
        Xtr = transform_scaler(Xtr_raw, means, stds)
        Xva = transform_scaler(Xva_raw, means, stds)

        best_lam, best_auc = lambda_grid[0], -1
        for lam in lambda_grid:
            w, b = fit_l1_logistic(Xtr, ytr, lam=lam)
            auc = auc_score(yva, predict_prob(Xva, w, b))
            if auc > best_auc:
                best_auc, best_lam = auc, lam
        return best_lam

    inner_k = min(4, max(2, n // 4))
    folds = kfold_indices(n, inner_k)
    best_lam, best_auc = lambda_grid[0], -1
    for lam in lambda_grid:
        val_pred = [0.5] * n
        for val_idx in folds:
            tr_idx = [i for i in range(n) if i not in set(val_idx)]
            if not tr_idx:
                continue
            Xtr_raw = [X_train_raw[i] for i in tr_idx]
            ytr = [y_train[i] for i in tr_idx]
            Xva_raw = [X_train_raw[i] for i in val_idx]

            means, stds = fit_scaler(Xtr_raw)
            Xtr = transform_scaler(Xtr_raw, means, stds)
            Xva = transform_scaler(Xva_raw, means, stds)

            w, b = fit_l1_logistic(Xtr, ytr, lam=lam)
            p = predict_prob(Xva, w, b)
            for idx, pp in zip(val_idx, p):
                val_pred[idx] = pp

        auc = auc_score(y_train, val_pred)
        if auc > best_auc:
            best_auc, best_lam = auc, lam
    return best_lam


def main():
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
    y = [1 if ph[s]['disease_group'] == cfg['PRIMARY_CASE'] else 0 for s in sample_ids]

    available = [g for g in feats if g in mat_symbol]
    if len(available) < 3:
        with open(RESULT_DIR / 'tables' / 'deg_primary_healthy_vs_npdr_pdr_dme.csv', encoding='utf-8') as f:
            top_ens = [r['gene'] for r in list(csv.DictReader(f))[:80]]
        rev = {k: v for k, v in mapping.items() if k in set(top_ens)}
        available = sorted(set(rev.values()))[:30]
        available = [g for g in available if g in mat_symbol]

    X_raw = []
    for s in sample_ids:
        idx = h[1:].index(s)
        X_raw.append([mat_symbol[g][idx] for g in available])

    outer_k = min(5, max(2, len(X_raw) // 4))
    folds = kfold_indices(len(X_raw), k=outer_k)
    lambda_grid = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2]

    oof = [0.5] * len(X_raw)
    chosen = []
    for test_idx in folds:
        train_idx = [i for i in range(len(X_raw)) if i not in set(test_idx)]
        Xtr_raw = [X_raw[i] for i in train_idx]
        ytr = [y[i] for i in train_idx]
        Xte_raw = [X_raw[i] for i in test_idx]

        best_lam = choose_lambda_inner_cv(Xtr_raw, ytr, lambda_grid)
        chosen.append(best_lam)

        means, stds = fit_scaler(Xtr_raw)
        Xtr = transform_scaler(Xtr_raw, means, stds)
        Xte = transform_scaler(Xte_raw, means, stds)

        w, b = fit_l1_logistic(Xtr, ytr, lam=best_lam)
        pte = predict_prob(Xte, w, b)
        for idx, pp in zip(test_idx, pte):
            oof[idx] = pp

    best_lam = sorted(chosen)[len(chosen) // 2]
    full_means, full_stds = fit_scaler(X_raw)
    X_full = transform_scaler(X_raw, full_means, full_stds)
    w_all, b_all = fit_l1_logistic(X_full, y, lam=best_lam)

    selected = [(g, c) for g, c in zip(available, w_all) if abs(c) > 1e-8]
    if len(selected) == 0:
        selected = sorted(zip(available, w_all), key=lambda x: -abs(x[1]))[:10]

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

    scores = []
    for i, s in enumerate(sample_ids):
        linear = b_all + sum(w_all[j] * X_full[i][j] for j in range(len(w_all)))
        scores.append({'sample_id': s, 'label': ph[s]['disease_group'], 'signature_score': linear, 'oof_probability': oof[i]})

    with open(RESULT_DIR / 'tables' / 'signature_scores.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=scores[0].keys())
        w.writeheader()
        w.writerows(scores)

    roc_rows = []
    for g in [x[0] for x in selected]:
        j = available.index(g)
        probs = [1.0 / (1.0 + math.exp(-X_full[i][j])) for i in range(len(X_full))]
        auc = auc_score(y, probs)
        lo, hi = bootstrap_auc_ci(y, probs)
        roc_rows.append({'item': g, 'auc': auc, 'ci95_low': lo, 'ci95_high': hi})

    with open(RESULT_DIR / 'tables' / 'roc_individual_genes.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['item', 'auc', 'ci95_low', 'ci95_high'])
        w.writeheader()
        w.writerows(sorted(roc_rows, key=lambda x: -x['auc']))

    auc = auc_score(y, oof)
    lo, hi = bootstrap_auc_ci(y, oof)
    with open(RESULT_DIR / 'tables' / 'roc_combined_signature.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['item', 'auc', 'ci95_low', 'ci95_high', 'lambda'])
        w.writeheader()
        w.writerow({'item': 'combined_signature_oof', 'auc': auc, 'ci95_low': lo, 'ci95_high': hi, 'lambda': best_lam})

    log_message('07_lasso_signature', f'features={len(available)} selected={len(selected)} auc={auc:.3f} lambda={best_lam}')


if __name__ == '__main__':
    main()
