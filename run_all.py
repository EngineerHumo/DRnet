import argparse, subprocess, json
from pathlib import Path
import runpy
from pipeline_utils import ensure_dirs
cfg=runpy.run_path('00_config.py')
RESULT_DIR=cfg['RESULT_DIR']
STEPS=['00a_validate_and_clean_gene_sets.py','01_parse_manifest.py','02_prepare_expression.py','03_qc_and_pca.py','04_differential_expression.py','05_inflammation_scoring.py','06_candidate_selection.py','07_lasso_signature.py','08_enrichment_analysis.py','09_immune_infiltration.py','10_make_figures.py','11_export_master_summary.py']

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--resume',action='store_true')
    args=ap.parse_args()
    ensure_dirs(RESULT_DIR/'logs')
    statef=RESULT_DIR/'logs'/'run_state.json'
    state={'completed':[]} if not statef.exists() else json.loads(statef.read_text())
    for step in STEPS:
        if args.resume and step in state['completed']:
            continue
        subprocess.run(['python',step],check=True)
        state['completed'].append(step)
        statef.write_text(json.dumps(state,indent=2),encoding='utf-8')
    print('Pipeline finished.')

if __name__=='__main__': main()
