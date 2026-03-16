import csv, os, math, json, random
from pathlib import Path
from collections import defaultdict, Counter

GROUPS = ["healthy control", "diabetic", "NPPR", "NPDR/PDR + DME"]
SEVERITY_MAP = {"healthy control":0,"diabetic":1,"NPPR":2,"NPDR/PDR + DME":3}

def ensure_dirs(*dirs):
    for d in dirs:
        Path(d).mkdir(parents=True, exist_ok=True)

def log_message(step, msg, logs_dir="results/logs"):
    ensure_dirs(logs_dir)
    with open(Path(logs_dir)/f"{step}.log","a",encoding="utf-8") as f:
        f.write(msg.rstrip()+"\n")

def read_gmt(path):
    sets={}
    with open(path,encoding='utf-8',errors='ignore') as f:
        for ln,line in enumerate(f,1):
            line=line.strip()
            if not line: continue
            parts=line.split('\t')
            if len(parts)<3:
                raise ValueError(f"Invalid GMT line {ln} in {path}")
            sets[parts[0]]=[g.strip() for g in parts[2:] if g.strip()]
    return sets

def write_gmt(path, sets):
    with open(path,'w',encoding='utf-8',newline='') as f:
        w=csv.writer(f,delimiter='\t')
        for k,v in sets.items():
            w.writerow([k,'na',*v])

def parse_series_matrix(path):
    rows={}
    with open(path,encoding='utf-8',errors='ignore') as f:
        for line in f:
            if not line.startswith('!Sample_'): continue
            rec=next(csv.reader([line],delimiter='\t'))
            key=rec[0].replace('!Sample_','')
            vals=[x.strip().strip('"') for x in rec[1:]]
            rows.setdefault(key,[]).append(vals)
    titles=rows['title'][0]
    data=[{'sample_id':s} for s in titles]
    for k,vv in rows.items():
        if k=='title':
            continue
        if k=='characteristics_ch1':
            for vals in vv:
                for i,v in enumerate(vals):
                    if ': ' in v:
                        kk,vv2=v.split(': ',1)
                        data[i][kk.strip()]=vv2.strip()
            continue
        vals=vv[0]
        for i,v in enumerate(vals):
            data[i][k]=v
    return data

def read_table(path):
    with open(path,encoding='utf-8',newline='') as f:
        r=csv.reader(f)
        header=next(r)
        rows=[row for row in r if row]
    return header,rows

def write_csv(path, rows, fieldnames):
    with open(path,'w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fieldnames)
        w.writeheader(); w.writerows(rows)

def bh_adjust(pvals):
    n=len(pvals)
    idx=sorted(range(n), key=lambda i:pvals[i])
    out=[1.0]*n
    prev=1.0
    for rank,i in enumerate(reversed(idx),1):
        k=n-rank+1
        val=min(prev,pvals[i]*n/k)
        out[i]=val; prev=val
    return out

def mann_whitney_u(x,y):
    vals=[(v,0) for v in x]+[(v,1) for v in y]
    vals.sort(key=lambda z:z[0])
    ranks=[0]*len(vals)
    i=0
    while i<len(vals):
        j=i
        while j<len(vals) and vals[j][0]==vals[i][0]: j+=1
        r=(i+j+1)/2
        for k in range(i,j): ranks[k]=r
        i=j
    rx=sum(r for r,(v,g) in zip(ranks,vals) if g==0)
    nx,ny=len(x),len(y)
    u=rx-nx*(nx+1)/2
    mu=nx*ny/2
    sigma=math.sqrt(nx*ny*(nx+ny+1)/12) if nx*ny>0 else 1
    z=(u-mu)/(sigma+1e-9)
    # normal approx two-sided
    p=2*(1-0.5*(1+math.erf(abs(z)/math.sqrt(2))))
    return u,p

def spearman(x,y):
    def rank(v):
        s=sorted((val,i) for i,val in enumerate(v))
        r=[0]*len(v);i=0
        while i<len(v):
            j=i
            while j<len(v) and s[j][0]==s[i][0]: j+=1
            rr=(i+j+1)/2
            for k in range(i,j): r[s[k][1]]=rr
            i=j
        return r
    rx,ry=rank(x),rank(y)
    mx,my=sum(rx)/len(rx),sum(ry)/len(ry)
    num=sum((a-mx)*(b-my) for a,b in zip(rx,ry))
    dx=math.sqrt(sum((a-mx)**2 for a in rx)); dy=math.sqrt(sum((b-my)**2 for b in ry))
    rho=num/(dx*dy+1e-12)
    n=len(x)
    t=abs(rho)*math.sqrt((n-2)/(1-rho*rho+1e-12)) if n>2 else 0
    p=2*(1-0.5*(1+math.erf(t/math.sqrt(2))))
    return rho,p

def ssgsea_score(expr_dict, geneset):
    genes=[g for g in geneset if g in expr_dict]
    if not genes:
        return 0.0
    return sum(expr_dict[g] for g in genes)/len(genes)
