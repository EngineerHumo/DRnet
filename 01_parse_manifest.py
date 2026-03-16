import csv, runpy
from collections import Counter
from pathlib import Path
from pipeline_utils import ensure_dirs, parse_series_matrix, write_csv, log_message
cfg=runpy.run_path('00_config.py')
RAW_DIR,PROC_DIR,RESULT_DIR=cfg['RAW_DIR'],cfg['PROC_DIR'],cfg['RESULT_DIR']
GROUPS=cfg['GROUPS'];PRIMARY_CTRL=cfg['PRIMARY_CTRL'];PRIMARY_CASE=cfg['PRIMARY_CASE']

def main():
    ensure_dirs(PROC_DIR, RESULT_DIR/'tables', RESULT_DIR/'logs', RAW_DIR)
    src=Path('GSE160306_series_matrix.txt'); dst=RAW_DIR/'GSE160306_series_matrix.txt'
    if src.exists() and not dst.exists(): dst.write_text(src.read_text(encoding='utf-8',errors='ignore'),encoding='utf-8')
    data=parse_series_matrix(dst)
    out=[]
    for r in data:
        raw_group=r.get('disease_group_dme', r.get('disease_group','')).strip()
        map_group={'Control':'healthy control','Diabetic':'diabetic','NPDR':'NPPR','NPDR/PDR + DME':'NPDR/PDR + DME'}
        g=map_group.get(raw_group,raw_group)
        region='macula' if 'Macula' in r.get('tissue','') else 'periphery'
        row={'sample_id':r['sample_id'],'geo_accession':r.get('geo_accession',''),'disease_group':g,'region':region,'severity_code':cfg['SEVERITY_MAP'].get(g,''),'disease_group_detailed':r.get('disease_group_detailed','')}
        out.append(row)
    write_csv(PROC_DIR/'manifest_all_samples.csv',out,list(out[0].keys()))
    mac=[r for r in out if r['region']=='macula' and r['disease_group'] in GROUPS]
    write_csv(PROC_DIR/'manifest_macula_4groups.csv',mac,list(mac[0].keys()))
    primary=[r for r in mac if r['disease_group'] in [PRIMARY_CTRL,PRIMARY_CASE]]
    for r in primary:r['primary_label']='case' if r['disease_group']==PRIMARY_CASE else 'control'
    write_csv(PROC_DIR/'manifest_primary_binary.csv',primary,list(primary[0].keys()))
    cnt=Counter((r['disease_group'],r['region']) for r in out)
    rows=[{'disease_group':k[0],'region':k[1],'n':v} for k,v in sorted(cnt.items())]
    write_csv(PROC_DIR/'sample_count_by_group_region.csv',rows,['disease_group','region','n'])
    log_message('01_parse_manifest',f'all={len(out)} macula={len(mac)} primary={len(primary)}')
if __name__=='__main__': main()
