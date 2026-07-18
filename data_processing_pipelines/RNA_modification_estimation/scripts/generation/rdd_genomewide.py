#!/usr/bin/env python3
"""
Genome-wide RNA-DNA difference (RDD) and substitution analysis from ONT
per-position pileups (the *_variation_tsv.gz files).

Reproduces the three-chromosome analysis (substitution matrix, RDD lower-bound,
per-modification directionality) across ALL chromosomes.

USAGE
-----
  # 1. Substitution matrix + RDD lower-bound (needs only the pileups):
  python3 rdd_genomewide.py matrix   PILEUP1.gz PILEUP2.gz ... 
  python3 rdd_genomewide.py rdd       PILEUP1.gz PILEUP2.gz ...

  # 2. Per-modification directionality (needs the bedRMod too):
  python3 rdd_genomewide.py permod  BEDRMOD.bed  PILEUP1.gz PILEUP2.gz ...

Pileups may be passed as a glob, e.g.:
  python3 rdd_genomewide.py rdd  /path/to/merged_native_noSupp_chr*_variation_tsv.gz

PARAMETERS (edit constants below to match the manuscript):
  MINCOV   = 20      minimum per-site coverage
  SNV_VAF  = 0.35    alt-fraction at/above which a site is treated as a genomic variant (excluded from RDDs)
  THRESHES = [...]   RDD level thresholds to scan
  Column indices assume the 23-column variation format:
    chrom,pos,ref,reads_all,reads_pp,matches,matches_pp,mismatches,mismatches_pp,
    deletions,deletions_pp,insertions,insertions_pp,A,A_pp,C,C_pp,T,T_pp,G,G_pp,N,N_pp
"""
import gzip, json, sys, glob
from collections import defaultdict, Counter

# ---------- parameters (match these to the manuscript) ----------
MINCOV   = 20
SNV_VAF  = 0.35
THRESHES = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
CONSENSUS = 17488          # genome-wide consensus, for reference printout
# ----------------------------------------------------------------

BASES = ['A','C','T','G']
COLIDX = {'A':13,'C':15,'T':17,'G':19}     # 0-based column of each base count
COMP  = {'A':'T','C':'G','G':'C','T':'A'}
DISP  = {'A':'A','C':'C','T':'U','G':'G'}  # display T as U
TRANSITIONS = {('A','G'),('G','A'),('C','T'),('T','C')}

def open_pileup(fn):
    f = gzip.open(fn,'rt')
    f.readline()  # skip header
    return f

def iter_sites(fn):
    """Yield (chrom,pos1,ref,counts_dict,total) for covered ACGT positions."""
    with open_pileup(fn) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 21: continue
            ref = p[2]
            if ref not in BASES: continue
            if int(p[3]) < MINCOV: continue
            cnt = {b:int(p[COLIDX[b]]) for b in BASES}
            tot = sum(cnt.values())
            if tot < MINCOV: continue
            yield p[0], int(p[1]), ref, cnt, tot

# ---------------- MODE 1: substitution matrix ----------------
def run_matrix(files):
    mat = {r:{b:0 for b in BASES} for r in BASES}
    refpos = {r:0 for r in BASES}
    for fn in files:
        sys.stderr.write(f"[matrix] {fn}\n")
        for chrom,pos,ref,cnt,tot in iter_sites(fn):
            refpos[ref]+=1
            for b in BASES: mat[ref][b]+=cnt[b]
    json.dump({'mat':mat,'refpos':refpos}, open('genomewide_matrix.json','w'))
    # report
    print("\n=== Substitution directionality (% of mismatches per ref base) ===")
    for r in BASES:
        mism={b:mat[r][b] for b in BASES if b!=r}; t=sum(mism.values())
        top=max(mism,key=mism.get)
        print(f"  {DISP[r]}: dominant {DISP[r]}>{DISP[top]} = {mism[top]/t*100:.0f}%  "
              + " ".join(f"{DISP[r]}>{DISP[b]} {mism[b]/t*100:.0f}%" for b in BASES if b!=r))
    print("\nSaved genomewide_matrix.json")

# ---------------- MODE 2: RDD lower-bound ----------------
def run_rdd(files):
    counts = {t:0 for t in THRESHES}
    dir_at = {t:Counter() for t in THRESHES}
    covered = 0
    for fn in files:
        sys.stderr.write(f"[rdd] {fn}\n")
        for chrom,pos,ref,cnt,tot in iter_sites(fn):
            covered+=1
            frac=(tot-cnt[ref])/tot
            if frac>=SNV_VAF: continue          # exclude genomic variants
            alt=max((b for b in BASES if b!=ref),key=lambda b:cnt[b])
            if (ref,alt) not in TRANSITIONS: continue   # transitions only
            for t in THRESHES:
                if frac>=t:
                    counts[t]+=1
                    dir_at[t][f"{DISP[ref]}>{DISP[alt]}"]+=1
    out={'thresh':{str(t):counts[t] for t in THRESHES},
         'covered':covered,
         'dir':{str(t):dict(dir_at[t]) for t in THRESHES}}
    json.dump(out, open('genomewide_rdd.json','w'))
    print("\n=== Transition RDDs genome-wide (excl variants >=%.0f%%, cov>=%d) ===" % (SNV_VAF*100, MINCOV))
    for t in THRESHES:
        print(f"  >={t:.0%}: {counts[t]:,}")
    print(f"\n  consensus reference: {CONSENSUS:,}")
    print(f"  covered positions: {covered:,}")
    print("\nSaved genomewide_rdd.json")

# ---------------- MODE 3: per-modification directionality ----------------
def run_permod(bedfile, files):
    # load mod calls keyed (chrom,pos1)->(mod,strand)
    modcall=defaultdict(dict)
    with open(bedfile) as f:
        for line in f:
            if line.startswith('#'): continue
            p=line.rstrip('\n').split('\t')
            modcall[p[0]][int(p[1])+1]=(p[3], p[5])   # bed start 0-based -> 1-based
    permod=defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    pos_count=Counter()
    for fn in files:
        sys.stderr.write(f"[permod] {fn}\n")
        with open_pileup(fn) as f:
            for line in f:
                p=line.rstrip('\n').split('\t')
                if len(p)<21: continue
                chrom=p[0]; 
                if chrom not in modcall: continue
                pos=int(p[1])
                if pos not in modcall[chrom]: continue
                ref=p[2]
                if ref not in BASES: continue
                if int(p[3])<MINCOV: continue
                cnt={b:int(p[COLIDX[b]]) for b in BASES}
                tot=sum(cnt.values())
                if tot<MINCOV: continue
                mod,strand=modcall[chrom][pos]
                if strand=='-':
                    ref=COMP[ref]; cnt={COMP[b]:cnt[b] for b in BASES}
                pos_count[mod]+=1
                for b in BASES: permod[mod][ref][b]+=cnt[b]
    out={m:{r:dict(permod[m][r]) for r in permod[m]} for m in permod}
    json.dump({'permod':out,'pos':dict(pos_count)}, open('genomewide_permod.json','w'))
    print("\n=== Per-modification directionality (RNA orientation) ===")
    modbase={'m6A':'A','I':'A','Am':'A','m5C':'C','Cm':'C','Y':'T','Um':'T','Gm':'G'}
    for m in sorted(pos_count,key=lambda k:-pos_count[k]):
        rb=modbase.get(m)
        if rb is None or rb not in permod[m]: continue
        row=permod[m][rb]; tot=sum(row.get(b,0) for b in BASES)
        if tot==0: continue
        mism={b:row.get(b,0) for b in BASES if b!=rb}; mt=sum(mism.values())
        top=max(mism,key=mism.get)
        print(f"  {m:>4} (n={pos_count[m]:>5}): {DISP[rb]}>{DISP[top]} = {mism[top]/mt*100:.0f}% of mismatches")
    print("\nSaved genomewide_permod.json")

# ---------------- dispatch ----------------
def expand(args):
    out=[]
    for a in args: out += sorted(glob.glob(a)) if any(c in a for c in '*?[') else [a]
    return out

if __name__=='__main__':
    if len(sys.argv)<3:
        print(__doc__); sys.exit(1)
    mode=sys.argv[1]
    if mode=='matrix': run_matrix(expand(sys.argv[2:]))
    elif mode=='rdd':  run_rdd(expand(sys.argv[2:]))
    elif mode=='permod':
        run_permod(sys.argv[2], expand(sys.argv[3:]))
    else:
        print("unknown mode:",mode); print(__doc__); sys.exit(1)
