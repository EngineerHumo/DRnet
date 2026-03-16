import csv, math, runpy
from collections import defaultdict
from pathlib import Path
from pipeline_utils import ensure_dirs, read_table, log_message
cfg=runpy.run_path('00_config.py')
RAW_DIR,PROC_DIR,RESULT_DIR=cfg['RAW_DIR'],cfg['PROC_DIR'],cfg['RESULT_DIR']

def read_manifest(path):
    with open(path,encoding='utf-8') as f:return list(csv.DictReader(f))

def subset_matrix(src, sample_ids):
    header,rows=read_table(src)
    idx=[0]+[header.index(s) for s in sample_ids if s in header]
    new_header=[header[i] for i in idx]
    out=[]
    for r in rows:
        out.append([r[i] for i in idx])
    return new_header,out

def main():
    ensure_dirs(PROC_DIR, RESULT_DIR/'logs')
    for fn in ['GSE160306_human_retina_DR_totalRNA_counts.txt','GSE160306_human_retina_DR_totalRNA_normalized_cpm.txt']:
        src=Path(fn); dst=RAW_DIR/fn
        if src.exists() and not dst.exists(): dst.write_text(src.read_text(encoding='utf-8',errors='ignore'),encoding='utf-8')
    mac=read_manifest(PROC_DIR/'manifest_macula_4groups.csv'); primary=read_manifest(PROC_DIR/'manifest_primary_binary.csv')
    mac_ids=[r['sample_id'] for r in mac]; pri_ids=[r['sample_id'] for r in primary]
    ch,mx=subset_matrix(RAW_DIR/'GSE160306_human_retina_DR_totalRNA_counts.txt',mac_ids)
    _,px=subset_matrix(RAW_DIR/'GSE160306_human_retina_DR_totalRNA_counts.txt',pri_ids)
    cpm_h,cpm=subset_matrix(RAW_DIR/'GSE160306_human_retina_DR_totalRNA_normalized_cpm.txt',mac_ids)
    # filter
    keep=[]
    for r in mx:
        vals=[float(x) for x in r[1:]]
        if sum(v>=cfg['MIN_COUNT'] for v in vals)>=cfg['MIN_SAMPLES']: keep.append(r[0])
    mx=[r for r in mx if r[0] in keep]; px=[r for r in px if r[0] in keep]; cpm=[r for r in cpm if r[0] in keep]
    for arr,name in [(mx,'counts_macula_4groups.tsv'),(px,'counts_primary_binary.tsv')]:
        with open(PROC_DIR/name,'w',newline='',encoding='utf-8') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(ch if 'macula' in name else [ch[0]]+pri_ids);w.writerows(arr)
    with open(PROC_DIR/'log2cpm_macula_4groups.tsv','w',newline='',encoding='utf-8') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(cpm_h)
        for r in cpm:w.writerow([r[0]]+[f"{math.log2(float(x)+1):.6f}" for x in r[1:]])
    for nm,data in [('pheno_macula_4groups.csv',mac),('pheno_primary_binary.csv',primary)]:
        with open(PROC_DIR/nm,'w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=data[0].keys());w.writeheader();w.writerows(data)
    with open(PROC_DIR/'gene_filtering_summary.csv','w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['genes_before','genes_after']);w.writeheader();w.writerow({'genes_before':len(read_table(RAW_DIR/'GSE160306_human_retina_DR_totalRNA_counts.txt')[1]),'genes_after':len(keep)})
    log_message('02_prepare_expression',f'genes_after_filter={len(keep)}')
if __name__=='__main__': main()
