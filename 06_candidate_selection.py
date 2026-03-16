import csv, runpy
from pipeline_utils import ensure_dirs, read_gmt, log_message
cfg=runpy.run_path('00_config.py')
RAW_DIR,RESULT_DIR=cfg['RAW_DIR'],cfg['RESULT_DIR']

def main():
    ensure_dirs(RESULT_DIR/'tables', RESULT_DIR/'logs')
    hall=set(read_gmt(RAW_DIR/'hallmark_inflammatory_response.gmt')['HALLMARK_INFLAMMATORY_RESPONSE'])
    with open(RESULT_DIR/'tables'/'deg_primary_healthy_vs_npdr_pdr_dme.csv',encoding='utf-8') as f:
        deg=[r for r in csv.DictReader(f) if int(r['significant'])==1]
    core=[{'gene':r['gene']} for r in deg if r['gene'] in hall]
    with open(RESULT_DIR/'tables'/'severity_trend_all_genes.csv',encoding='utf-8') as f:
        trend={r['gene']:r for r in csv.DictReader(f)}
    prog=[{'gene':r['gene'],'rho':trend[r['gene']]['spearman_rho'],'padj':trend[r['gene']]['padj']} for r in core if r['gene'] in trend and int(trend[r['gene']]['trend_significant'])==1]
    for nm,data in [('inflammatory_core_genes.csv',core),('progressive_inflammatory_genes.csv',prog)]:
        with open(RESULT_DIR/'tables'/nm,'w',newline='',encoding='utf-8') as f:
            if data:
                w=csv.DictWriter(f,fieldnames=data[0].keys());w.writeheader();w.writerows(data)
            else:
                f.write('gene\n')
    with open(RESULT_DIR/'tables'/'candidate_gene_summary.csv','w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['metric','value']);w.writeheader();w.writerows([{'metric':'primary_deg', 'value':len(deg)},{'metric':'inflammatory_core_genes','value':len(core)},{'metric':'progressive_inflammatory_genes','value':len(prog)}])
    log_message('06_candidate_selection',f'core={len(core)} progressive={len(prog)}')
if __name__=='__main__': main()
