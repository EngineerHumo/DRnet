import csv, math, runpy
from collections import defaultdict
from pipeline_utils import ensure_dirs, mann_whitney_u, bh_adjust, spearman, log_message
cfg=runpy.run_path('00_config.py')
PROC_DIR,RESULT_DIR=cfg['PROC_DIR'],cfg['RESULT_DIR']

def load_counts(path):
    with open(path,encoding='utf-8') as f:
        r=csv.reader(f,delimiter='\t');h=next(r);d={row[0]:[float(x) for x in row[1:]] for row in r}
    return h[1:],d

def run_contrast(samples, groups, data, g1,g2,out):
    i1=[i for i,s in enumerate(samples) if groups[s]==g1]; i2=[i for i,s in enumerate(samples) if groups[s]==g2]
    rows=[]
    for gene,vals in data.items():
        x=[vals[i] for i in i1]; y=[vals[i] for i in i2]
        lfc=math.log2((sum(y)/len(y)+1)/(sum(x)/len(x)+1))
        _,p=mann_whitney_u(x,y)
        rows.append({'gene':gene,'log2FC':lfc,'pvalue':p})
    padj=bh_adjust([r['pvalue'] for r in rows])
    for r,a in zip(rows,padj): r['padj']=a; r['significant']=int((a<cfg['PADJ_THRESHOLD']) and abs(r['log2FC'])>cfg['LOG2FC_THRESHOLD'])
    with open(out,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def main():
    ensure_dirs(RESULT_DIR/'tables', RESULT_DIR/'logs')
    samples,data=load_counts(PROC_DIR/'counts_macula_4groups.tsv')
    with open(PROC_DIR/'pheno_macula_4groups.csv',encoding='utf-8') as f: ph={r['sample_id']:r for r in csv.DictReader(f)}
    groups={k:v['disease_group'] for k,v in ph.items()}
    run_contrast(samples,groups,data,'healthy control','NPDR/PDR + DME',RESULT_DIR/'tables'/'deg_primary_healthy_vs_npdr_pdr_dme.csv')
    run_contrast(samples,groups,data,'healthy control','diabetic',RESULT_DIR/'tables'/'deg_diabetic_vs_healthy.csv')
    run_contrast(samples,groups,data,'healthy control','NPDR',RESULT_DIR/'tables'/'deg_nppr_vs_healthy.csv')
    run_contrast(samples,groups,data,'diabetic','NPDR/PDR + DME',RESULT_DIR/'tables'/'deg_npdr_pdr_dme_vs_diabetic.csv')
    sev={k:int(v['severity_code']) for k,v in ph.items()}
    trend=[]
    for gene,vals in data.items():
        x=[sev[s] for s in samples]; rho,p=spearman(x,vals)
        trend.append({'gene':gene,'spearman_rho':rho,'pvalue':p})
    padj=bh_adjust([r['pvalue'] for r in trend])
    for r,a in zip(trend,padj): r['padj']=a; r['trend_significant']=int(a<0.05 and r['spearman_rho']>0)
    with open(RESULT_DIR/'tables'/'severity_trend_all_genes.csv','w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=trend[0].keys());w.writeheader();w.writerows(trend)
    log_message('04_differential_expression','done de + trend')
if __name__=='__main__': main()
