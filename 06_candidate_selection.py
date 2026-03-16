import csv, runpy
from pipeline_utils import ensure_dirs, read_gmt, load_ensembl_symbol_mapping, normalize_ensembl_id, log_message

cfg = runpy.run_path('00_config.py')
RAW_DIR, PROC_DIR, RESULT_DIR = cfg['RAW_DIR'], cfg['PROC_DIR'], cfg['RESULT_DIR']


def main():
    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'logs')
    hallmark = set(read_gmt(RAW_DIR / 'hallmark_inflammatory_response.gmt')['HALLMARK_INFLAMMATORY_RESPONSE'])

    with open(RESULT_DIR / 'tables' / 'deg_primary_healthy_vs_npdr_pdr_dme.csv', encoding='utf-8') as f:
        deg = [r for r in csv.DictReader(f) if int(r['significant']) == 1]

    mapping = load_ensembl_symbol_mapping(RAW_DIR, PROC_DIR, ensembl_ids=[r['gene'] for r in deg])
    deg_symbol = []
    for r in deg:
        eid = normalize_ensembl_id(r['gene'])
        sym = mapping.get(eid)
        if sym:
            deg_symbol.append({'ensembl_id': eid, 'gene_symbol': sym})

    core = [r for r in deg_symbol if r['gene_symbol'] in hallmark]

    with open(RESULT_DIR / 'tables' / 'severity_trend_all_genes.csv', encoding='utf-8') as f:
        trend = {normalize_ensembl_id(r['gene']): r for r in csv.DictReader(f)}

    prog = []
    for r in core:
        t = trend.get(r['ensembl_id'])
        if t and int(t['trend_significant']) == 1:
            prog.append({
                'ensembl_id': r['ensembl_id'],
                'gene_symbol': r['gene_symbol'],
                'rho': t['spearman_rho'],
                'padj': t['padj']
            })

    for nm, data, fields in [
        ('inflammatory_core_genes.csv', core, ['ensembl_id', 'gene_symbol']),
        ('progressive_inflammatory_genes.csv', prog, ['ensembl_id', 'gene_symbol', 'rho', 'padj'])
    ]:
        with open(RESULT_DIR / 'tables' / nm, 'w', newline='', encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(data)

    with open(RESULT_DIR / 'tables' / 'candidate_gene_summary.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['metric', 'value'])
        w.writeheader()
        w.writerows([
            {'metric': 'primary_deg', 'value': len(deg)},
            {'metric': 'primary_deg_mapped_to_symbol', 'value': len(deg_symbol)},
            {'metric': 'inflammatory_core_genes', 'value': len(core)},
            {'metric': 'progressive_inflammatory_genes', 'value': len(prog)}
        ])

    with open(RESULT_DIR / 'logs' / '06_candidate_selection_mapping_check.md', 'w', encoding='utf-8') as f:
        f.write(
            f"# Candidate mapping check\n\n"
            f"- significant primary DEGs (ensembl): {len(deg)}\n"
            f"- mapped DEGs with gene symbol: {len(deg_symbol)}\n"
            f"- hallmark inflammatory genes: {len(hallmark)}\n"
            f"- inflammatory_core_genes: {len(core)}\n"
            f"- progressive_inflammatory_genes: {len(prog)}\n"
        )

    log_message('06_candidate_selection', f'core={len(core)} progressive={len(prog)} mapped_deg={len(deg_symbol)}')


if __name__ == '__main__':
    main()
