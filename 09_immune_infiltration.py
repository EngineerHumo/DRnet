import csv, runpy
from pipeline_utils import ensure_dirs, read_gmt, ssgsea_score, mann_whitney_u, spearman, bh_adjust, log_message
cfg=runpy.run_path('00_config.py')
PROC_DIR,RESULT_DIR=cfg['PROC_DIR'],cfg['RESULT_DIR']

def main():
    ensure_dirs(RESULT_DIR/'tables', RESULT_DIR/'logs')
    gmt=read_gmt('data_processed/immune_28_signatures_charoentong2017_human_cleaned.gmt')
    with open(PROC_DIR/'log2cpm_macula_4groups.tsv',encoding='utf-8') as f:
        r=csv.reader(f,delimiter='\t');h=next(r); mat={row[0]:[float(x) for x in row[1:]] for row in r}
    with open(PROC_DIR/'pheno_macula_4groups.csv',encoding='utf-8') as f: ph={r['sample_id']:r for r in csv.DictReader(f)}
    scores=[]
    for i,s in enumerate(h[1:]):
        expr={g:v[i] for g,v in mat.items()}
        for cell,genes in gmt.items():
            scores.append({'sample_id':s,'group':ph[s]['disease_group'],'severity_code':ph[s]['severity_code'],'cell_type':cell,'score':ssgsea_score(expr,genes)})
    with open(RESULT_DIR/'tables'/'immune_ssgsea_scores.csv','w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=scores[0].keys());w.writeheader();w.writerows(scores)
    cells=sorted(gmt)
    comp=[];trend=[]
    for c in cells:
        x=[r['score'] for r in scores if r['cell_type']==c and r['group']==cfg['PRIMARY_CTRL']]
        y=[r['score'] for r in scores if r['cell_type']==c and r['group']==cfg['PRIMARY_CASE']]
        _,p=mann_whitney_u(x,y); comp.append({'cell_type':c,'pvalue':p})
        rr=[r for r in scores if r['cell_type']==c]; rho,p2=spearman([int(a['severity_code']) for a in rr],[a['score'] for a in rr]); trend.append({'cell_type':c,'rho':rho,'pvalue':p2})
    for arr in [comp,trend]:
        adj=bh_adjust([r['pvalue'] for r in arr])
        for r,a in zip(arr,adj):r['padj']=a
    for nm,arr in [('immune_primary_comparison.csv',comp),('immune_severity_trend.csv',trend)]:
        with open(RESULT_DIR/'tables'/nm,'w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=arr[0].keys());w.writeheader();w.writerows(arr)
    with open(RESULT_DIR/'tables'/'lasso_selected_genes.csv',encoding='utf-8') as f: genes=[r['gene'] for r in csv.DictReader(f)]
    sig=[r['cell_type'] for r in comp if r['padj']<0.05]
    cor=[]
    for g in genes:
        if g not in mat: continue
        gv=mat[g]
        for c in sig:
            cv=[r['score'] for r in scores if r['cell_type']==c]
            rho,p=spearman(gv,cv)
            cor.append({'gene':g,'cell_type':c,'rho':rho,'pvalue':p})
    with open(RESULT_DIR/'tables'/'gene_immune_correlations.csv','w',newline='',encoding='utf-8') as f:
        fields=cor[0].keys() if cor else ['gene','cell_type','rho','pvalue']
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(cor)
    log_message('09_immune_infiltration',f'cells={len(cells)}')
if __name__=='__main__': main()
