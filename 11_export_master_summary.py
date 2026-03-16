import csv, runpy
from pipeline_utils import ensure_dirs, log_message
cfg=runpy.run_path('00_config.py')
RESULT_DIR=cfg['RESULT_DIR']

def main():
    ensure_dirs(RESULT_DIR/'tables', RESULT_DIR/'logs')
    files=[p for p in (RESULT_DIR/'tables').glob('*.csv')]
    rows=[]
    for p in sorted(files):
        with open(p,encoding='utf-8') as f:
            n=sum(1 for _ in f)-1
        rows.append({'table':p.name,'rows':max(n,0)})
    with open(RESULT_DIR/'tables'/'master_analysis_summary.csv','w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['table','rows']);w.writeheader();w.writerows(rows)
    with open(RESULT_DIR/'tables'/'master_analysis_summary.md','w',encoding='utf-8') as f:
        f.write('# Master analysis summary\n\n')
        for r in rows: f.write(f"- {r['table']}: {r['rows']} rows\n")
    log_message('11_export_master_summary',f'tables={len(rows)}')
if __name__=='__main__': main()
