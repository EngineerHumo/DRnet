import csv
import math
import random
import runpy
from pipeline_utils import ensure_dirs, read_gmt, load_ensembl_symbol_mapping, normalize_ensembl_id, bh_adjust, log_message

cfg = runpy.run_path('00_config.py')
RAW_DIR, PROC_DIR, RESULT_DIR = cfg['RAW_DIR'], cfg['PROC_DIR'], cfg['RESULT_DIR']
SEED = cfg.get('RANDOM_SEED', 202501)


def enrichment_score(ranked_genes, ranked_scores, gene_set, p=1.0):
    s = set(gene_set)
    hits = [i for i, g in enumerate(ranked_genes) if g in s]
    nh = len(hits)
    if nh == 0 or nh == len(ranked_genes):
        return 0.0
    nr = sum(abs(ranked_scores[i]) ** p for i in hits) + 1e-12
    miss = 1.0 / (len(ranked_genes) - nh)
    running = 0.0
    best_pos, best_neg = -1e9, 1e9
    hit_set = set(hits)
    for i in range(len(ranked_genes)):
        if i in hit_set:
            running += (abs(ranked_scores[i]) ** p) / nr
        else:
            running -= miss
        best_pos = max(best_pos, running)
        best_neg = min(best_neg, running)
    return best_pos if abs(best_pos) >= abs(best_neg) else best_neg


def perm_null(ranked_genes, ranked_scores, gene_set, n_perm=80):
    rnd = random.Random(SEED)
    scores = ranked_scores[:]
    null = []
    for _ in range(n_perm):
        rnd.shuffle(scores)
        null.append(enrichment_score(ranked_genes, scores, gene_set))
    return null


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

    ranked = sorted(gene_score.items(), key=lambda x: x[1], reverse=True)
    ranked_genes = [g for g, _ in ranked][:6000]
    ranked_scores = [s for _, s in ranked][:6000]

    hall = read_gmt(RAW_DIR / 'hallmark_human_gene_symbols.gmt')
    rows = []
    for pth, gset in hall.items():
        es = enrichment_score(ranked_genes, ranked_scores, gset)
        null = perm_null(ranked_genes, ranked_scores, gset)
        if es >= 0:
            pos = [x for x in null if x >= 0] or [1e-12]
            nes = es / (sum(pos) / len(pos))
            pval = sum(1 for x in pos if x >= es) / len(pos)
        else:
            neg = [x for x in null if x < 0] or [-1e-12]
            nes = -es / (abs(sum(neg) / len(neg)) + 1e-12)
            nes = -nes
            pval = sum(1 for x in neg if x <= es) / len(neg)
        hits = len(set(gset).intersection(set(ranked_genes)))
        rows.append({'pathway': pth, 'ES': es, 'NES': nes, 'pvalue': pval, 'hit_genes': hits, 'geneset_size': len(set(gset))})

    adj = bh_adjust([r['pvalue'] for r in rows])
    for r, a in zip(rows, adj):
        r['padj'] = a
    rows.sort(key=lambda x: (x['padj'], -abs(x['NES'])))

    with open(RESULT_DIR / 'tables' / 'gsea_summary_matrix.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    with open(RESULT_DIR / 'tables' / 'pathway_recurrence_summary.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['pathway', 'hit_genes', 'padj', 'NES'])
        w.writeheader()
        w.writerows([{'pathway': r['pathway'], 'hit_genes': r['hit_genes'], 'padj': r['padj'], 'NES': r['NES']} for r in rows])

    # gene-centric extraction from real GSEA result: pathways containing each selected gene, with global NES/FDR
    gene_rows = []
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
        gene_rows.extend(out)

    with open(RESULT_DIR / 'tables' / 'ipa_input_selected_genes.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['gene_symbol'])
        w.writeheader()
        for g in genes:
            w.writerow({'gene_symbol': g})

    log_message('08_enrichment_analysis', f'ranked_genes={len(ranked_genes)} pathways={len(rows)} selected_genes={len(genes)}')


if __name__ == '__main__':
    main()
