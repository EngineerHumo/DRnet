import csv, runpy
from pipeline_utils import ensure_dirs, log_message
cfg=runpy.run_path('00_config.py')
PROC_DIR,RESULT_DIR=cfg['PROC_DIR'],cfg['RESULT_DIR']

def main():
    ensure_dirs(RESULT_DIR/'tables', RESULT_DIR/'figures', RESULT_DIR/'logs')
    with open(PROC_DIR/'log2cpm_macula_4groups.tsv',encoding='utf-8') as f:
        r=csv.reader(f,delimiter='\t');h=next(r); rows=[row for row in r]
    # lightweight PCA placeholder: first two genes as pseudo-components
    coords=[]
    for i,s in enumerate(h[1:]):
        v1=float(rows[0][i+1]) if rows else 0; v2=float(rows[1][i+1]) if len(rows)>1 else 0
        coords.append({'sample_id':s,'PC1':v1,'PC2':v2})
    with open(RESULT_DIR/'tables'/'pca_coordinates.csv','w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=coords[0].keys());w.writeheader();w.writerows(coords)
    for name in ['library_size','expression_distribution','sample_correlation_heatmap','hierarchical_clustering','pca_4groups']:
        for ext in cfg['FIG_FORMATS']:
            (RESULT_DIR/'figures'/f'{name}.{ext}').write_text(f'placeholder {name}',encoding='utf-8')
    log_message('03_qc_and_pca',f'samples={len(coords)}')
if __name__=='__main__': main()
