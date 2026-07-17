#%%
import pandas as pd
import polars as pl
import os
import numpy as np
import matplotlib.pyplot as plt
import pyranges as pr
import seaborn as sns
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests
from upsetplot import UpSet, from_contents
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from upsetplot import UpSet,from_contents,plot
from matplotlib import pyplot
import re
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.metrics import jaccard_score
from scipy.stats import norm, rankdata
from tqdm import tqdm
import polars as pl
import math
from statsmodels.stats.multitest import multipletests
import io
import pyarrow
from joblib import Parallel, delayed
import matplotlib.gridspec as gridspec
from itertools import combinations
import matplotlib.colors as mcolors
from scipy.stats import rankdata, kendalltau
from concurrent.futures import ThreadPoolExecutor, as_completed



#%%


_MOD_MAP: dict[str, str] = {
    # --- Unmodified bases ---
    "A": "A",
    "16335": "A",                 # CHEBI:16335 — adenosine
    "G": "G",
    "16750": "G",                 # CHEBI:16750 — guanosine
    "C": "C",
    "17562": "C",                 # CHEBI:17562 — cytidine
    "U": "U",
    "16704": "U",                 # CHEBI:16704 — uridine
    "T": "T",

    # --- 2'-O-methyl ---
    "Am": "Am",
    "69426": "Am",                # CHEBI:69426 — 2'-O-methyladenosine
    "Cm": "Cm",
    "19228": "Cm",                # CHEBI:19228 — 2'-O-methylcytidine
    "Gm": "Gm",
    "19229": "Gm",                # CHEBI:19229 — 2'-O-methylguanosine
    "Um": "Um",
    "19227": "Um",                # CHEBI:19227 — 2'-O-methyluridine
    "Ym": "Ym",
    "m5Cm": "m5Cm",
    "184012": "m5Cm",             # CHEBI:184012 — 5,2'-O-dimethylcytidine
    "hm5Cm": "hm5Cm",
    "99997": "hm5Cm",             # placeholder (no public ChEBI ID)
    "f5Cm": "f5Cm",
    "99996": "f5Cm",              # placeholder (no public ChEBI ID)

    # --- Inosine ---
    "I": "I",
    "Ino": "I",
    "17596": "I",                 # CHEBI:17596 — inosine
    "A->G": "I",
    "m1I": "m1I",
    "19065": "m1I",               # CHEBI:19065 — 1-methylinosine

    # --- Pseudouridine ---
    "Y": "Y",
    "psU": "Y",
    "Psi": "Y",
    "17802": "Y",                 # CHEBI:17802 — pseudouridine
    "m1acp3Y": "m1acp3Y",
    "m1ap3U": "m1acp3Y",

    # --- Methyl-A ---
    "m1A": "m1A",
    "m1A ": "m1A",
    "16020": "m1A",               # CHEBI:16020 — 1-methyladenosine
    "m6A": "m6A",
    "a": "m6A",
    "x6A": "m6A",
    "21891": "m6A",               # CHEBI:21891 — N6-methyladenosine
    "m6,6A": "m6,6A",
    "m66A": "m6,6A",
    "28284": "m6,6A",             # CHEBI:28284 — N6,N6-dimethyladenosine
    "m7A": "m7A",

    # --- Methyl-C ---
    "m3C": "m3C",
    "m4C": "m4C",
    "m5C": "m5C",
    "m": "m5C",
    "20607": "m5C",               # CHEBI:20607 — 5-methylcytidine

    # --- Methyl-G ---
    "m1G": "m1G",
    "19062": "m1G",               # CHEBI:19062 — 1-methylguanosine
    "m2G": "m2G",
    "19702": "m2G",               # CHEBI:19702 — N2-methylguanosine
    "m2,2G": "m2,2G",
    "m22G": "m2,2G",
    "19289": "m2,2G",             # CHEBI:19289 — N2,N2-dimethylguanosine
    "m2,2,7G": "m2,2,7G",
    "143283": "m2,2,7G",          # CHEBI:143283 — N2,N2,7-trimethylguanosine
    "m6G": "m6G",
    "m7G": "m7G",
    "20794": "m7G",               # CHEBI:20794 — 7-methylguanosine

    # --- Methyl-U / dihydrouridine ---
    "m3U": "m3U",
    "89487": "m3U",               # CHEBI:89487 — 3-methyluridine
    "m5U": "m5U",
    "45996": "m5U",               # CHEBI:45996 — ribothymidine (5-methyluridine)
    "D": "D",
    "23774": "D",                 # CHEBI:23774 — dihydrouridine

    # --- Acetyl / acyl / threonyl ---
    "ac4C": "ac4C",
    "70989": "ac4C",              # CHEBI:70989 — N4-acetylcytidine
    "ac7G": "ac7G",
    "t6A": "t6A",
    "21440": "t6A",               # CHEBI:21440 — N6-threonylcarbamoyladenosine
    "m6t6A": "m6t6A",
    "133071": "m6t6A",            # CHEBI:133071 — N6-methyl-N6-threonylcarbamoyladenosine

    # --- C5-substituted U family ---
    "cm5U": "cm5U",
    "75654": "cm5U",              # CHEBI:75654 — 5-carboxymethyluridine
    "mcm5U": "mcm5U",
    "20598": "mcm5U",             # CHEBI:20598 — 5-methoxycarbonylmethyluridine
    "mcmo5U": "mcmo5U",
    "27241": "mcmo5U",            # CHEBI:27241 — 5-methoxycarbonylmethoxyuridine
    "mchm5U": "mchm5U",
    "99998": "mchm5U",            # placeholder (no public ChEBI ID)
    "0": "mchm5U",                # legacy placeholder
    "ncm5U": "ncm5U",
    "62005": "ncm5U",             # CHEBI:62005 — 5-carbamoylmethyluridine
    "ncm5Um": "ncm5Um",
    "99994": "ncm5Um",            # placeholder (no public ChEBI ID)
    "ncm5s2U": "ncm5s2U",
    "99995": "ncm5s2U",           # placeholder (no public ChEBI ID)

    # --- Hydroxymethyl / formyl C ---
    "hm5C": "hm5C",
    "191041": "hm5C",             # CHEBI:191041 — 5-hydroxymethylcytidine
    "f5C": "f5C",
    "234279": "f5C",              # CHEBI:234279 — 5-formylcytidine

    # --- Isopentenyl / thio A ---
    "i6A": "i6A",
    "62881": "i6A",               # CHEBI:62881 — N6-(Δ2-isopentenyl)adenosine
    "ms2i6A": "ms2i6A",
    "62875": "ms2i6A",            # CHEBI:62875 — 2-methylthio-N6-(Δ2-isopentenyl)adenosine
    "io6A": "io6A",
    "71693": "io6A",              # CHEBI:71693 — N6-(cis-hydroxyisopentenyl)adenosine
    "ms2io6A": "ms2io6A",
    "62879": "ms2io6A",           # CHEBI:62879 — 2-methylthio-N6-(cis-hydroxyisopentenyl)adenosine

    # --- 3-amino-3-carboxypropyl ---
    "acp3U": "acp3U",
    "acp3D": "acp3D",
    "71588": "acp3D",             # CHEBI:71588 — 3-(3-amino-3-carboxypropyl)dihydrouridine

    # --- Queuosine / 2-thio ---
    "Q": "Q",
    "60193": "Q",                 # CHEBI:60193 — queuosine
    "s2U": "s2U",
    "60731": "s2U",               # CHEBI:60731 — 2-thiouridine

    # --- Ambiguous / generic placeholders ---
    "mA?": "mA?",
    "99990": "mA?",               # placeholder
    "mC?": "mC?",
    "99991": "mC?",               # placeholder
    "mG?": "mG?",
    "99992": "mG?",               # placeholder
    "mU?": "mU?",
    "99993": "mU?",               # placeholder
    "U*": "U*",

    # --- Non-canonical / experimental placeholders ---
    "xp3Cm": "xp3Cm",
    "xp4U": "xp4U",
    "xp6G": "xp6G",
    "xp7G": "xp7G",
    "m2xp7G": "m2xp7G",

    # --- Base substitutions (not modifications) ---
    "G->T": "G->T",
    "G->A": "G->A",
    "G->C": "G->C",
    "C->A": "C->A",
    "C->T": "C->T",
    "C->G": "C->G",
    "A->C": "A->C",
    "A->T": "A->T",
    "T->C": "T->C",
    "T->G": "T->G",
    "T->A": "T->A",
}


_MOD_MAP_STRING = (
    "a:m6A:A,69426:Am:A,17596:I:A,Ino:I:A,I:I:A,A->G:I:A,"
    "17802:Y:U,Y:Y:U,A:A:A,G:G:G,C:C:C,U:U:U,T:m5U:U,"
    "m3C:m3C:C,psU:Y:U,Um:Um:U,Gm:Gm:G,Am:Am:A,Cm:Cm:C,"
    "19227:Um:U,19229:Gm:G,19228:Cm:C,"
    "m5C:m5C:C,20607:m5C:C,"
    "m6A:m6A:A,x6A:m6A:A,21891:m6A:A,"
    "xp3Cm:xp3Cm:C,m2xp7G:m2xp7G:G,m6G:m6G:G,xp4U:xp4U:U,ac7G:ac7G:G,"
    "89487:m3U:U,m3U:m3U:U,m5U:m5U:U,m4C:m4C:C,"
    "28284:m6,6A:A,xp6G:xp6G:G,xp7G:xp7G:G,m7A:m7A:A,"
    "m1acp3Y:m1acp3Y:U,m1ap3U:m1acp3Y:U,"
    "70989:ac4C:C,ac4C:ac4C:C,"
    "20794:m7G:G,m7G:m7G:G,m66A:m6,6A:A,Ym:Ym:U,"
    "16020:m1A:A,m1A:m1A:A,m1A :m1A:A,m:m5C:C,"
    "143283:m2,2,7G:G,"
    "G->T::G,T->C::U,G->A::G,C->A::C,C->T::C,A->C::A,"
    "G->C::G,T->G::U,T->A::U,C->G::C,A->T::A,"
    "mA?:mA?:A,mC?:mC?:C,mG?:mG?:G,mU?:mU?:U,"
    # --- previously missing ---
    # unmodified CHEBI IDs
    "16335:A:A,16750:G:G,17562:C:C,16704:U:U,"
    # 2'-O-methyl C5 variants
    "m5Cm:m5Cm:C,184012:m5Cm:C,"
    "hm5Cm:hm5Cm:C,99997:hm5Cm:C,"
    "f5Cm:f5Cm:C,99996:f5Cm:C,"
    # 1-methylinosine
    "m1I:m1I:A,19065:m1I:A,"
    # pseudouridine alias
    "Psi:Y:U,"
    # m6,6A canonical name
    "m6,6A:m6,6A:A,"
    # methylguanosines
    "m1G:m1G:G,19062:m1G:G,"
    "m2G:m2G:G,19702:m2G:G,"
    "m2,2G:m2,2G:G,m22G:m2,2G:G,19289:m2,2G:G,"
    "m2,2,7G:m2,2,7G:G,"
    # ribothymidine / dihydrouridine
    "45996:m5U:U,"
    "D:D:U,23774:D:U,"
    # threonylcarbamoyl
    "t6A:t6A:A,21440:t6A:A,"
    "m6t6A:m6t6A:A,133071:m6t6A:A,"
    # C5-substituted U family
    "cm5U:cm5U:U,75654:cm5U:U,"
    "mcm5U:mcm5U:U,20598:mcm5U:U,"
    "mcmo5U:mcmo5U:U,27241:mcmo5U:U,"
    "mchm5U:mchm5U:U,99998:mchm5U:U,0:mchm5U:U,"
    "ncm5U:ncm5U:U,62005:ncm5U:U,"
    "ncm5Um:ncm5Um:U,99994:ncm5Um:U,"
    "ncm5s2U:ncm5s2U:U,99995:ncm5s2U:U,"
    # hydroxymethyl / formyl C
    "hm5C:hm5C:C,191041:hm5C:C,"
    "f5C:f5C:C,234279:f5C:C,"
    # isopentenyl A family
    "i6A:i6A:A,62881:i6A:A,"
    "ms2i6A:ms2i6A:A,62875:ms2i6A:A,"
    "io6A:io6A:A,71693:io6A:A,"
    "ms2io6A:ms2io6A:A,62879:ms2io6A:A,"
    # acp3 family
    "acp3U:acp3U:U,"
    "acp3D:acp3D:U,71588:acp3D:U,"
    # queuosine / 2-thiouridine
    "Q:Q:G,60193:Q:G,"
    "s2U:s2U:U,60731:s2U:U,"
    # placeholder numeric IDs
    "99990:mA?:A,99991:mC?:C,99992:mG?:G,99993:mU?:U,"
    # misc
    "U*:U*:U"
)



_CHROM_MAP: dict[str, str] = {
    # **{i:f"chr{i}" for i in range(1, 23)},
    **{str(i): f"chr{i}" for i in range(1, 23)},
    "X": "chrX",
    "Y": "chrY",
    "MT": "chrM",
    **{f"chr{i}": f"chr{i}" for i in range(1, 23)},
    "chrX": "chrX",
    "chrY": "chrY",
    "chrM": "chrM",
    "chrMT": "chrM",
    #

    # GRCh38
    "NC_000001.11": "chr1",
    "NC_000002.12": "chr2",
    "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",
    "NC_000005.10": "chr5",
    "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",
    "NC_000008.11": "chr8",
    "NC_000009.12": "chr9",
    "NC_000010.11": "chr10",
    "NC_000011.10": "chr11",
    "NC_000012.12": "chr12",
    "NC_000013.11": "chr13",
    "NC_000014.09": "chr14",
    "NC_000015.10": "chr15",
    "NC_000016.10": "chr16",
    "NC_000017.11": "chr17",
    "NC_000018.10": "chr18",
    "NC_000019.10": "chr19",
    "NC_000020.11": "chr20",
    "NC_000021.09": "chr21",
    "NC_000022.11": "chr22",
    "NC_000023.11": "chrX",
    "NC_000024.10": "chrY",
    "NC_012920.1": "chrM",
    # GRCh37
    "NC_000001.10": "chr1",
    "NC_000002.11": "chr2",
    "NC_000003.11": "chr3",
    "NC_000004.11": "chr4",
    "NC_000005.9": "chr5",
    "NC_000006.11": "chr6",
    "NC_000007.13": "chr7",
    "NC_000008.10": "chr8",
    "NC_000009.11": "chr9",
    "NC_000010.10": "chr10",
    "NC_000011.9": "chr11",
    "NC_000012.11": "chr12",
    "NC_000013.10": "chr13",
    "NC_000014.8": "chr14",
    "NC_000015.9": "chr15",
    "NC_000016.9": "chr16",
    "NC_000017.10": "chr17",
    "NC_000018.9": "chr18",
    "NC_000019.9": "chr19",
    "NC_000020.10": "chr20",
    "NC_000021.8": "chr21",
    "NC_000022.10": "chr22",
    "NC_000023.10": "chrX",
    "NC_000024.9": "chrY",
    "NC_001807.4": "chrM",
    # scaffolds
    # GRCh38 unplaced/unlocalized scaffolds
    "GL000008.2": "GL000008.2",
    "GL000009.2": "GL000009.2",
    "GL000194.1": "GL000194.1",
    "GL000195.1": "GL000195.1",
    "GL000205.2": "GL000205.2",
    "GL000208.1": "GL000208.1",
    "GL000213.1": "GL000213.1",
    "GL000214.1": "GL000214.1",
    "GL000216.2": "GL000216.2",
    "GL000218.1": "GL000218.1",
    "GL000219.1": "GL000219.1",
    "GL000220.1": "GL000220.1",
    "GL000221.1": "GL000221.1",
    "GL000224.1": "GL000224.1",
    "GL000225.1": "GL000225.1",
    "GL000226.1": "GL000226.1",
    "KI270302.1": "KI270302.1",
    "KI270303.1": "KI270303.1",
    "KI270304.1": "KI270304.1",
    "KI270305.1": "KI270305.1",
    "KI270310.1": "KI270310.1",
    "KI270311.1": "KI270311.1",
    "KI270312.1": "KI270312.1",
    "KI270315.1": "KI270315.1",
    "KI270316.1": "KI270316.1",
    "KI270317.1": "KI270317.1",
    "KI270320.1": "KI270320.1",
    "KI270322.1": "KI270322.1",
    "KI270329.1": "KI270329.1",
    "KI270330.1": "KI270330.1",
    "KI270333.1": "KI270333.1",
    "KI270334.1": "KI270334.1",
    "KI270335.1": "KI270335.1",
    "KI270336.1": "KI270336.1",
    "KI270337.1": "KI270337.1",
    "KI270338.1": "KI270338.1",
    "KI270340.1": "KI270340.1",
    "KI270362.1": "KI270362.1",
    "KI270363.1": "KI270363.1",
    "KI270364.1": "KI270364.1",
    "KI270366.1": "KI270366.1",
    "KI270371.1": "KI270371.1",
    "KI270372.1": "KI270372.1",
    "KI270373.1": "KI270373.1",
    "KI270374.1": "KI270374.1",
    "KI270375.1": "KI270375.1",
    "KI270376.1": "KI270376.1",
    "KI270378.1": "KI270378.1",
    "KI270379.1": "KI270379.1",
    "KI270381.1": "KI270381.1",
    "KI270382.1": "KI270382.1",
    "KI270383.1": "KI270383.1",
    "KI270384.1": "KI270384.1",
    "KI270385.1": "KI270385.1",
    "KI270386.1": "KI270386.1",
    "KI270387.1": "KI270387.1",
    "KI270388.1": "KI270388.1",
    "KI270389.1": "KI270389.1",
    "KI270390.1": "KI270390.1",
    "KI270391.1": "KI270391.1",
    "KI270392.1": "KI270392.1",
    "KI270393.1": "KI270393.1",
    "KI270394.1": "KI270394.1",
    "KI270395.1": "KI270395.1",
    "KI270396.1": "KI270396.1",
    "KI270411.1": "KI270411.1",
    "KI270412.1": "KI270412.1",
    "KI270414.1": "KI270414.1",
    "KI270417.1": "KI270417.1",
    "KI270418.1": "KI270418.1",
    "KI270419.1": "KI270419.1",
    "KI270420.1": "KI270420.1",
    "KI270422.1": "KI270422.1",
    "KI270423.1": "KI270423.1",
    "KI270424.1": "KI270424.1",
    "KI270425.1": "KI270425.1",
    "KI270429.1": "KI270429.1",
    "KI270435.1": "KI270435.1",
    "KI270438.1": "KI270438.1",
    "KI270442.1": "KI270442.1",
    "KI270448.1": "KI270448.1",
    "KI270465.1": "KI270465.1",
    "KI270466.1": "KI270466.1",
    "KI270467.1": "KI270467.1",
    "KI270468.1": "KI270468.1",
    "KI270507.1": "KI270507.1",
    "KI270508.1": "KI270508.1",
    "KI270509.1": "KI270509.1",
    "KI270510.1": "KI270510.1",
    "KI270511.1": "KI270511.1",
    "KI270512.1": "KI270512.1",
    "KI270515.1": "KI270515.1",
    "KI270516.1": "KI270516.1",
    "KI270517.1": "KI270517.1",
    "KI270518.1": "KI270518.1",
    "KI270519.1": "KI270519.1",
    "KI270521.1": "KI270521.1",
    "KI270522.1": "KI270522.1",
    "KI270528.1": "KI270528.1",
    "KI270529.1": "KI270529.1",
    "KI270530.1": "KI270530.1",
    "KI270538.1": "KI270538.1",
    "KI270539.1": "KI270539.1",
    "KI270544.1": "KI270544.1",
    "KI270548.1": "KI270548.1",
    "KI270579.1": "KI270579.1",
    "KI270580.1": "KI270580.1",
    "KI270581.1": "KI270581.1",
    "KI270582.1": "KI270582.1",
    "KI270583.1": "KI270583.1",
    "KI270584.1": "KI270584.1",
    "KI270587.1": "KI270587.1",
    "KI270588.1": "KI270588.1",
    "KI270589.1": "KI270589.1",
    "KI270590.1": "KI270590.1",
    "KI270591.1": "KI270591.1",
    "KI270593.1": "KI270593.1",
    "KI270706.1": "KI270706.1",
    "KI270707.1": "KI270707.1",
    "KI270708.1": "KI270708.1",
    "KI270709.1": "KI270709.1",
    "KI270710.1": "KI270710.1",
    "KI270711.1": "KI270711.1",
    "KI270712.1": "KI270712.1",
    "KI270713.1": "KI270713.1",
    "KI270714.1": "KI270714.1",
    "KI270715.1": "KI270715.1",
    "KI270716.1": "KI270716.1",
    "KI270717.1": "KI270717.1",
    "KI270718.1": "KI270718.1",
    "KI270719.1": "KI270719.1",
    "KI270720.1": "KI270720.1",
    "KI270721.1": "KI270721.1",
    "KI270722.1": "KI270722.1",
    "KI270723.1": "KI270723.1",
    "KI270724.1": "KI270724.1",
    "KI270725.1": "KI270725.1",
    "KI270726.1": "KI270726.1",
    "KI270727.1": "KI270727.1",
    "KI270728.1": "KI270728.1",
    "KI270729.1": "KI270729.1",
    "KI270730.1": "KI270730.1",
    "KI270731.1": "KI270731.1",
    "KI270732.1": "KI270732.1",
    "KI270733.1": "KI270733.1",
    "KI270734.1": "KI270734.1",
    "KI270735.1": "KI270735.1",
    "KI270736.1": "KI270736.1",
    "KI270737.1": "KI270737.1",
    "KI270738.1": "KI270738.1",
    "KI270739.1": "KI270739.1",
    "KI270740.1": "KI270740.1",
    "KI270741.1": "KI270741.1",
    "KI270742.1": "KI270742.1",
    "KI270743.1": "KI270743.1",
    "KI270744.1": "KI270744.1",
    "KI270745.1": "KI270745.1",
    "KI270746.1": "KI270746.1",
    "KI270747.1": "KI270747.1",
    "KI270748.1": "KI270748.1",
    "KI270749.1": "KI270749.1",
    "KI270750.1": "KI270750.1",
    "KI270751.1": "KI270751.1",
    "KI270752.1": "KI270752.1",
    "KI270753.1": "KI270753.1",
    "KI270754.1": "KI270754.1",
    "KI270755.1": "KI270755.1",
    "KI270756.1": "KI270756.1",
    "KI270757.1": "KI270757.1",
    "GL000251.2": "GL000251.2", #Not in hg38 combined ref
    "GL000252.2": "GL000252.2", #Not in hg38 combined ref
    "GL000253.2": "GL000253.2", #Not in hg38 combined ref
    "GL000254.2": "GL000254.2", #Not in hg38 combined ref
    "GL000255.2": "GL000255.2", #Not in hg38 combined ref
    "GL000256.2": "GL000256.2", #Not in hg38 combined ref
    "GL949747.2": "GL949747.2", #Not in hg38 combined ref
    "GL949748.2": "GL949748.2", #Not in hg38 combined ref
    "GL949749.2": "GL949749.2", #Not in hg38 combined ref
    "GL949750.2": "GL949750.2", #Not in hg38 combined ref
    "GL949751.2": "GL949751.2", #Not in hg38 combined ref
    "GL949753.2": "GL949753.2", #Not in hg38 combined ref
    "KI270776.1": "KI270776.1", #Not in hg38 combined ref
    "KI270816.1": "KI270816.1", #Not in hg38 combined ref
    "KI270853.1": "KI270853.1", #Not in hg38 combined ref
    "KI270879.1": "KI270879.1", #Not in hg38 combined ref
    "KI270908.1": "KI270908.1", #Not in hg38 combined ref
    "KN196479.1": "KN196479.1", #Not in hg38 combined ref
    "KN538360.1": "KN538360.1", #Not in hg38 combined ref
    "KN538361.1": "KN538361.1", #Not in hg38 combined ref
    "KQ090018.1": "KQ090018.1", #Not in hg38 combined ref
    "KV575244.1": "KV575244.1", #Not in hg38 combined ref
    "KZ208914.1": "KZ208914.1", #Not in hg38 combined ref
    "ML143343.1": "ML143343.1", #Not in hg38 combined ref
    "MU273333.1": "MU273333.1", #Not in hg38 combined ref
    "MU273376.1": "MU273376.1", #Not in hg38 combined ref
    "chr22_KI270733v1_random": "chr22_KI270733v1_random", #Not in hg38 combined ref
    "chr1_GL383518v1_alt": "chr1_GL383518v1_alt", #Not in hg38 combined ref
    "chr2_KI270767v1_alt": "chr2_KI270767v1_alt", #Not in hg38 combined ref
    "chr5_KI270791v1_alt": "chr5_KI270791v1_alt", #Not in hg38 combined ref
    "chr6_GL000250v2_alt": "chr6_GL000250v2_alt", #Not in hg38 combined ref
    "chr7_KI270803v1_alt": "chr7_KI270803v1_alt", #Not in hg38 combined ref
    "chr15_KI270850v1_alt": "chr15_KI270850v1_alt", #Not in hg38 combined ref
    "chr16_KI270853v1_alt": "chr16_KI270853v1_alt", #Not in hg38 combined ref
    "chr17_KI270857v1_alt": "chr17_KI270857v1_alt", #Not in hg38 combined ref
    "chr19_GL383574v1_alt": "chr19_GL383574v1_alt", #Not in hg38 combined ref
    "chr2_KI270893v1_alt": "chr2_KI270893v1_alt", #Not in hg38 combined ref
    "chr6_GL000251v2_alt": "chr6_GL000251v2_alt", #Not in hg38 combined ref
    "chr6_GL000252v2_alt": "chr6_GL000252v2_alt", #Not in hg38 combined ref
    "chr6_GL000253v2_alt": "chr6_GL000253v2_alt", #Not in hg38 combined ref
    "chr6_GL000254v2_alt": "chr6_GL000254v2_alt", #Not in hg38 combined ref
    "chr6_GL000255v2_alt": "chr6_GL000255v2_alt", #Not in hg38 combined ref
    "chr6_GL000256v2_alt": "chr6_GL000256v2_alt", #Not in hg38 combined ref
    "chrUn_GL000220v1": "chrUn_GL000220v1", #Not in hg38 combined ref
    "chr14_KI270847v1_alt": "chr14_KI270847v1_alt", #Not in hg38 combined ref
    "chr14_KI270846v1_alt": "chr14_KI270846v1_alt", #Not in hg38 combined ref
    "chr16_KI270855v1_alt": "chr16_KI270855v1_alt", #Not in hg38 combined ref


    # rRNA
    "hs_rRNA_18S": "hs_rRNA_18S",
    "hs_rRNA_28S": "hs_rRNA_28S",
    "hs_rRNA_5.8S": "hs_rRNA_5.8S",
    "hs_rRNA_5S": "hs_rRNA_5S",
    "18S": "hs_rRNA_18S",
    "28S": "hs_rRNA_28S",
    "5.8S": "hs_rRNA_5.8S",
    "5S": "hs_rRNA_5S",
    "hs_mt-rRNA_12s": "hs_mt-rRNA_12s",
    "hs_mt-rRNA_16s": "hs_mt-rRNA_16s",
    # tRNA
    'hs_tRNAAla_CGC': 'hs_tRNAAla_CGC',
    'hs_tRNAAla_TGC': 'hs_tRNAAla_TGC',
    'hs_tRNAAla_AGC': 'hs_tRNAAla_AGC',
    'hs_tRNAArg_ACG': 'hs_tRNAArg_ACG',
    'hs_tRNAArg_TCG': 'hs_tRNAArg_TCG',
    'hs_tRNAArg_CCG2': 'hs_tRNAArg_CCG2',
    'hs_tRNAArg_TCT2': 'hs_tRNAArg_TCT2',
    'hs_tRNAArg_CCT': 'hs_tRNAArg_CCT',
    'hs_tRNAArg_CCG1': 'hs_tRNAArg_CCG1',
    'hs_tRNAArg_TCT1': 'hs_tRNAArg_TCT1',
    'hs_tRNAAsn_GTT': 'hs_tRNAAsn_GTT',
    'hs_tRNAAsp_GTC': 'hs_tRNAAsp_GTC',
    'hs_tRNACys_GCA': 'hs_tRNACys_GCA',
    'hs_tRNAGln_CTG_TTG': 'hs_tRNAGln_CTG_TTG',
    'hs_tRNAGlu_CTC': 'hs_tRNAGlu_CTC',
    'hs_tRNAGlu_TTC': 'hs_tRNAGlu_TTC',
    'hs_tRNAGly_GCC': 'hs_tRNAGly_GCC',
    'hs_tRNAGly_TCC': 'hs_tRNAGly_TCC',
    'hs_tRNAGly_CCC': 'hs_tRNAGly_CCC',
    'hs_tRNAHis_GTG': 'hs_tRNAHis_GTG',
    'hs_tRNAIle_AAT': 'hs_tRNAIle_AAT',
    'hs_tRNAIle_TAT': 'hs_tRNAIle_TAT',
    'hs_tRNAIle_GAT': 'hs_tRNAIle_GAT',
    'hs_tRNALeu_AAG': 'hs_tRNALeu_AAG',
    'hs_tRNALeu_CAG': 'hs_tRNALeu_CAG',
    'hs_tRNALeu_CAA': 'hs_tRNALeu_CAA',
    'hs_tRNALeu_TAA': 'hs_tRNALeu_TAA',
    'hs_tRNALys_CTT': 'hs_tRNALys_CTT',
    'hs_tRNALys_TTT': 'hs_tRNALys_TTT',
    'hs_tRNAMet_CAT': 'hs_tRNAMet_CAT',
    'hs_tRNAPhe_GAA': 'hs_tRNAPhe_GAA',
    'hs_tRNAPro_AGG_CGG_TGG': 'hs_tRNAPro_AGG_CGG_TGG',
    'hs_tRNASeC_TCA': 'hs_tRNASeC_TCA',
    'hs_tRNASer_AGA': 'hs_tRNASer_AGA',
    'hs_tRNASer_CGA': 'hs_tRNASer_CGA',
    'hs_tRNASer_GCT': 'hs_tRNASer_GCT',
    'hs_tRNAThr_AGT': 'hs_tRNAThr_AGT',
    'hs_tRNAThr_CGT1': 'hs_tRNAThr_CGT1',
    'hs_tRNAThr_TGT': 'hs_tRNAThr_TGT',
    'hs_tRNAThr_CGT2': 'hs_tRNAThr_CGT2',
    'hs_tRNATrp_CCA': 'hs_tRNATrp_CCA',
    'hs_tRNATyr_GTA2': 'hs_tRNATyr_GTA2',
    'hs_tRNATyr_ATA': 'hs_tRNATyr_ATA',
    'hs_tRNATyr_GTA1': 'hs_tRNATyr_GTA1',
    'hs_tRNAVal_AAC_CAC': 'hs_tRNAVal_AAC_CAC',
    'hs_tRNAVal_TAC': 'hs_tRNAVal_TAC',
    'hs_tRNAArg_CCG': 'hs_tRNAArg_CCG1',
    'hs_tRNAArg_TCT': 'hs_tRNAArg_TCT1',
    'hs-tRNAArg-TCG': 'hs_tRNAArg_TCG',
    'hs_tRNATyr_GTA' : 'hs_tRNATyr_GTA1',
    'hs_tRNAiMet_CAT': 'hs_tRNAiMet_CAT',
    'hs_tRNALeu_TAG': 'hs_tRNALeu_TAG', #Not in hg38 combined ref
    'hs_tRNASer_TGA': 'hs_tRNASer_TGA', #Not in hg38 combined ref
    'hs_tRNAThr_CGT': 'hs_tRNAThr_CGT1', 
    'hs_mttRNAAla_TGC': 'hs_mttRNAAla_TGC',
    'hs_mttRNAArg_TCG': 'hs_mttRNAArg_TCG',
    'hs_mttRNAAsn_GTT': 'hs_mttRNAAsn_GTT',
    'hs_mttRNAAsp_GTC': 'hs_mttRNAAsp_GTC',
    'hs_mttRNACys_GCA': 'hs_mttRNACys_GCA',
    'hs_mttRNAGln_TTG': 'hs_mttRNAGln_TTG',
    'hs_mttRNAGlu_TTC': 'hs_mttRNAGlu_TTC',
    'hs_mttRNAGly_TCC': 'hs_mttRNAGly_TCC',
    'hs_mttRNAHis_GTG': 'hs_mttRNAHis_GTG',
    'hs_mttRNAIle_GAT': 'hs_mttRNAIle_GAT',
    'hs_mttRNALeu_TAA': 'hs_mttRNALeu_TAA',
    'hs_mttRNALeu_TAG': 'hs_mttRNALeu_TAG',
    'hs_mttRNALys_TTT': 'hs_mttRNALys_TTT',
    'hs_mttRNAMet_CAT': 'hs_mttRNAMet_CAT',
    'hs_mttRNAPhe_GAA': 'hs_mttRNAPhe_GAA',
    'hs_mttRNAPro_TGG': 'hs_mttRNAPro_TGG',
    'hs_mttRNASer_GCT': 'hs_mttRNASer_GCT',
    'hs_mttRNASer_TGA': 'hs_mttRNASer_TGA',
    'hs_mttRNAThr_TGT': 'hs_mttRNAThr_TGT',
    'hs_mttRNATrp_TCA': 'hs_mttRNATrp_TCA',
    'hs_mttRNATyr_GTA': 'hs_mttRNATyr_GTA',
    'hs_mttRNAVal_TAC': 'hs_mttRNAVal_TAC',
    'hs_tRNAVal_AAC_CAC_G34=I_introduced': 'hs_tRNAVal_AAC_CAC',
    'hs_tRNASer_AGA_G34=I_introduced': 'hs_tRNASer_AGA',
    'hs_tRNAPro_AGG_CGG_TGG_G34=I_introduced': 'hs_tRNAPro_AGG_CGG_TGG',
    'hs_tRNALeu_AAG_G34=I_introduced': 'hs_tRNALeu_AAG',
    'hs_tRNAArg_ACG_G34=I_introduced': 'hs_tRNAArg_ACG',
    'hs_tRNAAla_AGC_G34=I_introduced': 'hs_tRNAAla_AGC',
    'hs_tRNAVal_AAC_CAC G34=I_introduced': 'hs_tRNAVal_AAC_CAC',
    'hs_tRNASer_AGA G34=I_introduced': 'hs_tRNASer_AGA',
    'hs_tRNAPro_AGG_CGG_TGG G34=I_introduced': 'hs_tRNAPro_AGG_CGG_TGG',
    'hs_tRNALeu_AAG G34=I_introduced': 'hs_tRNALeu_AAG',
    'hs_tRNAArg_ACG G34=I_introduced': 'hs_tRNAArg_ACG',
    'hs_tRNAAla_AGC G34=I_introduced': 'hs_tRNAAla_AGC',
    
    'oligo3' : 'oligo3',
    'oligo5' : 'oligo5',
    'oligo5_oligo3' : 'oligo5_oligo3'
    }


def shift_coordinates(df: pl.DataFrame, offset: int = 0) -> pl.DataFrame:
    """
    Shifts the start and end coordinate columns of a BED-style DataFrame by a given offset.

    Parameters:
        df (pl.DataFrame): Input DataFrame where the second column is the start coordinate
                           and the third column is the end coordinate.
        offset (int): Integer value to add to both start and end coordinates. Defaults to 0.

    Returns:
        pl.DataFrame: DataFrame with updated start and end coordinate columns.
    """
    start, end = df.columns[1], df.columns[2]
    return df.with_columns(
        (pl.col(start) + offset).alias(start),
        (pl.col(end) + offset).alias(end),
    )


# ── Scalar helpers (kept for any non-DataFrame call sites) ───────────────────
def modbasecode_to_mod_polars(code: str) -> str:
    """
    Converts a modification base code to its standardised modification name.

    Parameters:
        code (str): The modification base code to look up.

    Returns:
        str: The corresponding standardised modification name from _MOD_MAP.
    """
    return _MOD_MAP[str(code)]


def chromname_renaming_polars(code: str) -> str:
    """
    Converts a chromosome name to its UCSC-style equivalent (e.g. '1' -> 'chr1').

    Parameters:
        code (str): The chromosome name to look up.

    Returns:
        str: The mapped UCSC-style chromosome name, or the original code if not found in _CHROM_MAP.
    """
    return _CHROM_MAP.get(code, code)


# ── Polars-native DataFrame functions ────────────────────────────────────────


def align_modification_names_polars(df: pl.DataFrame) -> pl.DataFrame:
    """
    Maps the 'name' column to standardised modification names using _MOD_MAP.

    Unmapped values are kept as-is. The 'name' column is cast to String before mapping.

    Parameters:
        df (pl.DataFrame): Input DataFrame containing a 'name' column with modification codes.

    Returns:
        pl.DataFrame: DataFrame with the 'name' column replaced by standardised modification names.
    """
    return df.with_columns(
        pl.col("name")
        .cast(pl.String)
        .replace(_MOD_MAP, default=pl.col("name"))
        .alias("name")
    )


def align_chromosome_names_polars(df: pl.DataFrame) -> pl.DataFrame:
    """
    Maps the 'chrom' column to UCSC-style chromosome names (e.g. '1' -> 'chr1').

    Values not present in _CHROM_MAP are passed through unchanged.

    Parameters:
        df (pl.DataFrame): Input DataFrame containing a 'chrom' column.

    Returns:
        pl.DataFrame: DataFrame with the 'chrom' column mapped to UCSC-style names.
    """
    for i in df["chrom"]:
        try:
            _CHROM_MAP[i]
        except:
            print(i)
    return df.with_columns(
        pl.col("chrom")
        .cast(pl.String)
        .replace_strict(
            old=list(_CHROM_MAP.keys()),
            new=list(_CHROM_MAP.values()),
            return_dtype=pl.String
        )
        .alias("chrom")
    )


def load_bedrmod_polars(
    path: str = "", modtype: str = "", with_comment: bool = True
) -> tuple[pl.DataFrame, dict]:
    """
    Loads a bedRMod file into a Polars DataFrame.

    Parses comment lines (prefixed with '#') into a key-value dictionary, assigns
    standardised column names based on column count, normalises the frequency column
    to a percentage scale, aligns modification and chromosome names, and filters out
    rows with zero frequency or empty modification names.

    Parameters:
        path (str): Path to the bedRMod file.
        modtype (str): If provided, overrides all values in the 'name' column with
                       this modification type string.
        with_comment (bool): Whether to skip comment lines when reading the file.
                             Defaults to True.

    Returns:
        tuple[pl.DataFrame, dict]:
            - DataFrame containing the first 11 standard bedRMod columns.
            - Dictionary of key-value pairs parsed from the file's comment lines.
    """
    COLUMN_MAP = {
        11: [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "coverage",
            "frequency",
        ],
        12: [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "coverage",
            "frequency",
            "unknown",
        ],
        13: [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "coverage",
            "frequency",
            "rmsk",
            "db",
        ],
        14: [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "coverage",
            "frequency",
            "unique_mapping",
            "frag_start",
            "frag_end",
        ],
        ##chrom	chromStart	chromEnd	name	strand	score	thickStart	thickEnd	itemRgb	coverage	frequency	frequency_std	scores_std	overlap_perc	overlap
        15: [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "strand",
            "score",
            "thickStart",
            "thickEnd",	
            "itemRgb",
            "coverage",
            "frequency",
            "frequency_std",
            "scores_std",
            "overlap_perc",
            "overlap"
        ],
        17:["chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "coverage",
            "frequency",
            "nRep",
            "repName", 
            "repScore",
            "repCov",  
            "repFreq method"
        ],
        18: [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "coverage",
            "frequency",
            "count_modified",
            "count_canonical",
            "count_other_mod",
            "count_delete",
            "count_fail",
            "count_diff",
            "count_nocall",
        ],
        19: [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "coverage",
            "frequency",
            "n_mod",
            "count_canonical",
            "count_other_mod",
            "count_delete",
            "count_fail",
            "count_diff",
            "count_nocall",
            "score_nmod_thresholds",
        ],
        23: [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "coverage",
            "frequency",
            "n_mod",
            "count_canonical",
            "count_other_mod",
            "count_delete",
            "count_fail",
            "count_diff",
            "count_nocall",
            "score_nmod_thresholds",
            "std_coverage",
            "std_n_mod",
            "std_score_nmod_thresholds",
            "std_frequency",
        ],
        24: [
            "chrom", 
            "chromStart", 
            "chromEnd", 
            "name", 
            "score", 
            "strand", 
            "thickStart", 
            "thickEnd", 
            "itemRgb", 
            "coverage", 
            "frequency", 
            "n_mod", 
            "count_canonical", 
            "count_other_mod", 
            "count_delete", 
            "count_fail", 
            "count_diff",
            "count_nocall", 
            "score_nmod_thresholds", 
            "std_coverage", 
            "std_n_mod", 
            "std_score_nmod_thresholds", 
            "std_frequency", 
            "type"
        ],
        28: [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "coverage",
            "frequency",
            "n_mod",
            "count_canonical",
            "count_other_mod",
            "count_delete",
            "count_fail",
            "count_diff",
            "count_nocall",
            "score_nmod_thresholds",
            "std_coverage",
            "std_n_mod",
            "std_score_nmod_thresholds",
            "std_frequency",
            "pvalue_shuffle_axis_1",
            "pvalue_adj_shuffle_axis_1",
            "pvalue_shuffle_axis_0",
            "pvalue_ranks",
            "pvalue_shuffle_axis_1_IVT",
        ],
    }

    comment_dict: dict = {}
    raw_lines: list[str] = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#"):
                stripped = line.strip()
                key_value = stripped[1:].split("=", 1)
                if len(key_value) == 2:
                    key, value = key_value
                    comment_dict[key.strip()] = value.strip()
            else:
                raw_lines.append(line)

    df = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        comment_prefix="#" if with_comment else None,
        infer_schema_length=10000,
    )
    n_cols = len(df.columns)
    if n_cols in COLUMN_MAP:
        df = df.rename(dict(zip(df.columns, COLUMN_MAP[n_cols])))
    else:
        print(f"Warning: unexpected column count ({n_cols}); columns left unnamed.")

    if df["frequency"].max() <= 1.0:
        df = df.with_columns((pl.col("frequency") * 100).alias("frequency"))


    if modtype:
        df = df.with_columns(pl.lit(modtype).alias("name"))

    df = align_modification_names_polars(df)
    df = align_chromosome_names_polars(df)

    try:
        df = df.filter(pl.col("frequency") > 0)
    except Exception:
        print("Frequency is not numeric")

    df = df.filter(pl.col("name") != "")

    return df, comment_dict


def concatenate_bed_comments(
    bed_comments: list[dict], rna_types: list[str], technology: str
) -> dict:
    """
    Merges per-file bedRMod comment dictionaries into a single combined comment dictionary.

    Each metadata field (assembly, annotation_source, basecalling, etc.) is aggregated
    across all files using the pattern '<rna_type>:<value>;' so that the origin of each
    entry is traceable in the output file header.

    Parameters:
        bed_comments (list[dict]): List of comment dictionaries, one per input bedRMod file,
                                   as returned by load_bedrmod_polars.
        rna_types (list[str]): List of RNA type labels corresponding to each file
                               (e.g. ['mRNA', 'tRNA']). Must be the same length as bed_comments.
        technology (str): Sequencing platform/technology string written to the
                          'sequencing_plattform' field of the output header.

    Returns:
        dict: A merged comment dictionary suitable for passing to write_bedrmod_polars.
    """
    concatenated_dict = {}
    concatenated_dict["fileformat"] = "bedRModv2"
    concatenated_dict["modification_type"] = "".join(
        f"{str(r_type)}\t" for r_type in rna_types
    )
    concatenated_dict["modification_names"] = _MOD_MAP_STRING
    concatenated_dict["sequencing_plattform"] = technology
    assembly_string = ""
    annotation_source_string = ""
    basecalling_string = ""
    bioinformatics_workflow_string = ""
    experiment_string = ""
    external_source_string = ""
    for comment, rna_type in zip(bed_comments, rna_types):
        assembly_string += f"{rna_type}:{comment.get('assembly', '')};"
        annotation_source_string += (
            f"{rna_type}:{comment.get('annotation_source', '')};"
        )
        basecalling_string += f"{rna_type}:{comment.get('basecalling', '')};"
        bioinformatics_workflow_string += (
            f"{rna_type}:{comment.get('bioinformatics_workflow', '')};"
        )
        experiment_string += f"{rna_type}:[{comment.get('experiment', '')}];"
        external_source_string += f"{rna_type}:{comment.get('external_source', '')};"
    concatenated_dict["assembly"] = assembly_string
    concatenated_dict["annotation_source"] = annotation_source_string
    concatenated_dict["basecalling"] = basecalling_string
    concatenated_dict["bioinformatics_workflow"] = bioinformatics_workflow_string
    concatenated_dict["experiment"] = experiment_string
    concatenated_dict["external_source"] = external_source_string
    return concatenated_dict


def write_bedrmod_polars(
    input_df: pl.DataFrame, input_comment_dict: dict, output_filename: str
) -> None:
    """
    Writes a Polars DataFrame to a bedRMod-formatted file.

    Comment key-value pairs are written first as '#key=value' header lines,
    followed by a tab-separated column header line prefixed with '#', and then
    the data rows without a header.

    Parameters:
        input_df (pl.DataFrame): DataFrame to write. Column names are used as the
                                 header line in the output file.
        input_comment_dict (dict): Ordered dictionary of metadata to write as comment
                                   lines at the top of the file.
        output_filename (str): Destination file path for the output bedRMod file.

    Returns:
        None
    """
    with open(output_filename, "w") as file:
        for key, value in input_comment_dict.items():
            file.write(f"#{key}={value}\n")
        column_header = "".join(f"{str(col)}\t" for col in input_df.columns)
        file.write(f"#{column_header}\n")
        input_df.write_csv(file, separator="\t", include_header=False)







def calculate_bedrmod_score_weighted_sum(
    data: pl.DataFrame,
    data2: pl.DataFrame = None,
    n_permutations: int = 1000,
    min_coverage: int = 1,
    null_sample_size: int = 500,   # <-- rows sampled per permutation
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Compute a weighted-sum BEDRMod score for RNA modification sites, with
    empirical p-values derived from a subsampled permutation null distribution
    and Benjamini–Hochberg multiple-testing correction.

    The pipeline consists of six steps:

    1. **Coverage filtering** – retain only positions with coverage ≥
       ``min_coverage``.
    2. **Min-max standardisation** – coverage and n_mod are log₂-transformed
       before scaling; frequency is divided by 100; score_nmod_thresholds is
       scaled directly. All four features are mapped to [0, 1].
    3. **Weighted sum** – the four standardised features are combined with equal
       weights (0.25 each) to produce an observed score per position.
    4. **Permutation null distribution** – for each permutation, ``null_sample_size``
       rows are drawn without replacement and their feature matrix is shuffled
       row-wise before scoring, building a pooled null distribution of size
       ``n_permutations × null_sample_size``.
    5. **Empirical p-values** – computed via binary search on the sorted null
       pool, with a +1 pseudocount to avoid p = 0.
    6. **BH correction** – adjusted p-values are computed at α = 0.05, and both
       raw and adjusted BEDRMod scores are appended to the result.

    Parameters
    ----------
    data : pl.DataFrame
        Input Polars DataFrame. Must contain at minimum the columns
        ``"frequency"`` (0–100 float), ``"coverage"`` (int), and
        ``"score_nmod_thresholds"`` (float). Positions with ``frequency == 0``
        are dropped before any further processing.
    n_permutations : int, default 1000
        Number of permutation iterations used to build the null distribution.
        Higher values improve p-value resolution at the cost of runtime.
    min_coverage : int, default 1
        Minimum read coverage required to retain a position. Positions below
        this threshold are excluded before scoring.
    null_sample_size : int, default 500
        Number of rows sampled (without replacement) per permutation. Capped
        internally at the number of available positions. Reducing this value
        decreases memory usage (null pool ≈ ``n_permutations × null_sample_size``
        floats) while introducing additional variance in the null distribution.
    verbose : bool, default True
        If ``True``, prints the number of positions retained after coverage
        filtering and the total null pool size.

    Returns
    -------
    pd.DataFrame
        The filtered DataFrame extended with the following columns:

        - ``std_coverage``, ``std_n_mod``, ``std_score_nmod_thresholds``,
          ``std_frequency`` – standardised feature values.
        - ``pvalue`` – empirical p-value from the permutation null distribution.
        - ``pvalue_adj`` – BH-adjusted p-value.
        - ``score`` – BEDRMod score derived from the raw p-value.
        - ``score_adj`` – BEDRMod score derived from the adjusted p-value.

    Raises
    ------
    ValueError
        If ``"frequency"`` or ``"coverage"`` columns are absent from ``data``.
    ValueError
        If any computed ``n_mod`` value is negative or exceeds ``coverage``.
    ValueError
        If no positions remain after applying the ``min_coverage`` filter.

    Warns
    -----
    UserWarning
        If fewer than 100 positions pass the coverage filter, indicating that
        p-value estimates may be unreliable.

    Notes
    -----
    - ``n_mod`` is derived internally as ``floor(frequency × coverage / 100)``.
    - The permutation null shuffles rows within each sampled sub-matrix (i.e.
      feature vectors are permuted across positions), breaking the association
      between genomic loci and their feature profiles.
    - Empirical p-values use the formula ``(count_geq + 1) / (n_null + 1)``
      (Phipson & Smyth, 2010) to ensure they are strictly bounded above zero.
    - The random number generator is seeded at 42 for reproducibility.
    - Time complexity of the p-value step is O(N log M), where N is the number
      of positions and M = ``n_permutations × null_sample_size``.
    """
    def _minmax_expr(col: str) -> pl.Expr:
        if col == "coverage" or col == "n_mod":
            c = (pl.col(col) + 1).log(base=2)
        else:
            c = pl.col(col)
        return ((c - c.min()) / (c.max() - c.min())).alias(f"std_{col}")
    # --- Convert to Polars ---
    df = data
    if not data2.is_empty():
        df2 = data2
        df2 = df2.filter(pl.col("frequency") > 0).with_columns(
        ((pl.col("frequency") * pl.col("coverage")) // 100)
        .cast(pl.Int64)
        .alias("n_mod")
        )
        df2 = df2.with_columns([
        _minmax_expr("coverage").alias("std_coverage"),
        _minmax_expr("n_mod").alias("std_n_mod"),
        _minmax_expr("score_nmod_thresholds").alias("std_score_nmod_thresholds"),
        (pl.col("frequency") / 100).alias("std_frequency"),
        ])
        df2 = df2.with_columns((pl.col("std_coverage") * 0.3333 + 
                        pl.col("std_frequency") * 0.3333 + 
                        pl.col("std_score_nmod_thresholds") * 0.3333).alias("score"))

    # --- Input validation ---
    required = {"frequency", "coverage"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Data must contain columns: {', '.join(missing)}")

    df = df.filter(pl.col("frequency") > 0).with_columns(
        ((pl.col("frequency") * pl.col("coverage")) // 100)
        .cast(pl.Int64)
        .alias("n_mod")
    )

    

    if (df["n_mod"] < 0).any() or (df["coverage"] < 0).any():
        raise ValueError("n_mod and coverage cannot be negative")
    if (df["n_mod"] > df["coverage"]).any():
        raise ValueError("n_mod cannot exceed coverage")

    # --- Step 1: Filter coverage ---
    n_orig = len(df)
    df = df.filter(pl.col("coverage") >= min_coverage)

    if df.is_empty():
        raise ValueError(f"No positions pass min_coverage threshold of {min_coverage}")
    if len(df) < 100:
        warnings.warn("Less than 100 positions after filtering – p-values may be unreliable")
    if verbose:
        print(f"Retained {len(df)} of {n_orig} positions after min_coverage filter")

    # --- Step 2: Min-max standardisation ---


    df = df.with_columns([
        _minmax_expr("coverage").alias("std_coverage"),
        _minmax_expr("n_mod").alias("std_n_mod"),
        _minmax_expr("score_nmod_thresholds").alias("std_score_nmod_thresholds"),
        (pl.col("frequency") / 100).alias("std_frequency"),
    ])

    # df = df.with_columns((pl.col("std_coverage") * 0.25 + 
    #                       pl.col("std_frequency") * 0.25 + 
    #                       pl.col("std_n_mod") * 0.25 + 
    #                       pl.col("std_score_nmod_thresholds") * 0.25).alias("weighted_sum"))
    df = df.with_columns((pl.col("std_coverage") * 0.3333 + 
                          pl.col("std_frequency") * 0.3333 + 
                          pl.col("std_score_nmod_thresholds") * 0.3333).alias("score"))
    

    # --- Step 3: Weighted sum ---
    # std_cols = ["std_coverage", "std_n_mod", "std_score_nmod_thresholds", "std_frequency"]
    # weights  = np.array([0.25, 0.25, 0.25, 0.25])
    std_cols = ["std_coverage", "std_score_nmod_thresholds", "std_frequency"]
    weights  = np.array([0.3333, 0.3333, 0.3333])

    mat      = df.select(std_cols).to_numpy() 
    mat_IVT = df2.select(std_cols).to_numpy()  # (N, 4)
    observed = mat @ weights                     # (N,)

    # --- Step 4.1: Subsampled permutation null distribution ---
    N = len(df)
    k = min(null_sample_size, N)
    rng = np.random.default_rng(42)
 
    null_pool = np.empty(n_permutations * k, dtype=np.float64)
 
    for i in tqdm(range(n_permutations)):
        # Sample k rows from the full feature matrix
        idx = rng.choice(N, size=k, replace=False)
        null_mat = mat[idx].copy()
        for col in range(null_mat.shape[1]):
            rng.shuffle(null_mat[:, col])
 
        null_pool[i * k : (i + 1) * k] = null_mat @ weights

    # --- Step 5.1: Empirical p-values via searchsorted  ---
    # Runtime:(O(N log M))
    null_sorted = np.sort(null_pool)              # ascending
    n_null      = len(null_sorted)
    count_geq   = n_null - np.searchsorted(null_sorted, observed, side="left")
    pvalues     = (count_geq + 1) / (n_null + 1)  # +1 pseudocount, avoids p=0

    if verbose:
        print(
            f"Null pool size: {n_null:,}  "
            f"(null_sample_size={k}, n_permutations={n_permutations})"
        )

    # --- Step 6.1: BH correction + bedrmod score ---
    _, pvalue_adj, _, _ = multipletests(pvalues, method="fdr_bh", alpha=0.05)

    df = df.with_columns([
        pl.Series("pvalue_shuffle_axis_1",     pvalues),
        pl.Series("pvalue_adj_shuffle_axis_1", pvalue_adj)
    ])

    # --- Step 4.2: Subsampled permutation null distribution ---
    N = len(df)
    # Cap sample size so it never exceeds available rows
    k = min(null_sample_size, N)
    rng = np.random.default_rng(42)

    # Pool all null scores into a single 1-D array
    # Memory: (n_permutations * k,) floats
    null_pool_2 = np.empty(n_permutations * k, dtype=np.float64)

    for i in tqdm(range(n_permutations)):
        idx = rng.choice(N, size=k, replace=False)
        # Permute only the sampled rows before scoring
        null_pool_2[i * k : (i + 1) * k] = rng.permuted(mat[idx], axis=0) @ weights

    # --- Step 5.2: Empirical p-values via searchsorted  ---
    # Runtime:(O(N log M))
    null_sorted_2 = np.sort(null_pool_2)              # ascending
    n_null      = len(null_sorted_2)
    count_geq   = n_null - np.searchsorted(null_sorted_2, observed, side="left")
    pvalues     = (count_geq + 1) / (n_null + 1)  # +1 pseudocount, avoids p=0

    df = df.with_columns([
        pl.Series("pvalue_shuffle_axis_0",     pvalues)
    ])

    # --- Step 4.3 & 5.3: Pvalue from ranks ---
    ranks   = rankdata(-observed, method="average")
    pvalues = ranks / (N + 1)

    df = df.with_columns([
        pl.Series("pvalue_ranks",     pvalues)
    ])
    print(df.head())
    print(df2.head())

    if not df2.is_empty():
        # ---Step 4.4 P-vlaue when considering IVT dataset as null ---
        N = len(df2)
        k = min(null_sample_size, N)
        rng = np.random.default_rng(42)
    
        null_pool_3 = np.empty(n_permutations * k, dtype=np.float64)
    
        for i in tqdm(range(n_permutations)):
            # Sample k rows from the full feature matrix
            idx = rng.choice(N, size=k, replace=False)
            null_mat = mat_IVT[idx].copy()
            for col in range(null_mat.shape[1]):
                rng.shuffle(null_mat[:, col])
    
            null_pool_3[i * k : (i + 1) * k] = null_mat @ weights

        # --- Step 5.4: Empirical p-values via searchsorted  ---
        # Runtime:(O(N log M))
        null_sorted_3 = np.sort(null_pool_3)              # ascending
        n_null      = len(null_sorted_3)
        count_geq   = n_null - np.searchsorted(null_sorted_3, observed, side="left")
        pvalues     = (count_geq + 1) / (n_null + 1)  # +1 pseudocount, avoids p=0

        df = df.with_columns([
            pl.Series("pvalue_shuffle_axis_1_IVT",     pvalues)
        ])
        # ---Step 4.5 P-vlaue when considering IVT dataset as null ---
        mat_IVT = df2.select(std_cols).to_numpy() 
        N = len(df2)
        k = min(null_sample_size, N)
        rng = np.random.default_rng(42)
    
        null_pool_4 = np.empty(n_permutations * k, dtype=np.float64)
    
        for i in tqdm(range(n_permutations)):
            # Sample k rows from the full feature matrix
            idx = rng.choice(N, size=k, replace=False)
            null_mat = mat_IVT[idx].copy()
            # Permute only the sampled rows before scoring
            null_pool_4[i * k : (i + 1) * k] = rng.permuted(null_mat, axis=0) @ weights
        
        # --- Step 5.5: Empirical p-values via searchsorted  ---
        # Runtime:(O(N log M))
        null_sorted_4 = np.sort(null_pool_4)              # ascending
        n_null      = len(null_sorted_4)
        count_geq   = n_null - np.searchsorted(null_sorted_4, observed, side="left")
        pvalues     = (count_geq + 1) / (n_null + 1)  # +1 pseudocount, avoids p=0
        df = df.with_columns([
            pl.Series("pvalue_shuffle_axis_0_IVT",     pvalues)
        ])
    return df, df2     #, null_sorted, null_sorted_2, null_sorted_3



#%%

names_of_dfs= []
list_of_native_dfs = []
list_of_native_dfs_comments = []
list_of_IVT_dfs = []
list_of_IVT_dfs_comments = []
list_of_indices = []

#Load all files from threshold 0.85-0.99
for confidence_index,i in enumerate(range(85,100,1)):
    string_i = f"0.{i}"
    native_df_i, native_df_comment = load_bedrmod_polars(f"{string_i}.bed")

    native_df_i.columns = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","count_modified","count_canonical", "count_other_mod", "count_delete", "count_fail", "count_diff", "count_nocall"]
    native_df_i = native_df_i.filter(~pl.col("chrom").str.contains("rRNA"))
    native_df_i = native_df_i.filter(~pl.col("chrom").str.contains("oligo"))
    native_df_i = native_df_i.filter(pl.col("frequency") > 0)
    native_df_i = native_df_i.with_columns(pl.lit(confidence_index+1).alias("score_nmod_thresholds"))
    print(native_df_i.shape)
    
    ivt_df_i, ivt_df_comment = load_bedrmod_polars(f"{string_i}.bed")
    ivt_df_i.columns = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","count_modified","count_canonical", "count_other_mod", "count_delete", "count_fail", "count_diff", "count_nocall"]
    ivt_df_i = ivt_df_i.filter(~pl.col("chrom").str.contains("rRNA"))
    ivt_df_i = ivt_df_i.filter(~pl.col("chrom").str.contains("oligo"))
    ivt_df_i = ivt_df_i.filter(pl.col("frequency") > 0)
    ivt_df_i = ivt_df_i.with_columns(pl.lit(confidence_index+1).alias("score_nmod_thresholds"))
    print(ivt_df_i.shape)
    
    names_of_dfs.append(string_i)
    list_of_native_dfs.append(native_df_i)
    list_of_native_dfs_comments.append(native_df_comment)
    list_of_IVT_dfs.append(ivt_df_i)
    list_of_IVT_dfs_comments.append(ivt_df_comment)
    list_of_indices.append(confidence_index + 1)
    


# %%
#Load 0.85 modification confidence score table as initial dataframe
final_native_df = pl.DataFrame(list_of_native_dfs[0])
final_ivt_df = pl.DataFrame(list_of_IVT_dfs[0])


# Update maximal model confidence score so that the threshold boundary where modification occur is reach (0.85,0.99 in 0.01 stepsize)
def update_score(df_source: pl.DataFrame, df_target: pl.DataFrame) -> pl.DataFrame:
    keys = ["chrom", "chromStart", "name", "score", "strand"]
    return (
        df_source
        .join(
            df_target.select(keys + ["score_nmod_thresholds"]),
            on=keys,
            how="left",
            suffix="_new",
        )
        .with_columns(
            pl.coalesce("score_nmod_thresholds_new", "score_nmod_thresholds")
            .alias("score_nmod_thresholds")
        )
        .drop("score_nmod_thresholds_new")
    )




# Run update score iteratively to update the maximal model confidence boundary
for native_df,ivt_df in zip(list_of_native_dfs, list_of_IVT_dfs):
    final_native_df = update_score(final_native_df,native_df)
    final_ivt_df = update_score(final_ivt_df,ivt_df)

# Calculate the weghted sum score used with downstream filtering --> 0.33 * (min_max(log2(coverage)) + (frequency / 100) + min_max(score_nmod_thresholds))
final_native_df, final_ivt_df = calculate_bedrmod_score_weighted_sum(data = final_native_df,data2 = final_ivt_df)    
    
"""
The following section should only run with tRNA data, 
since the tRNA sequencing was conducted with a 24nt long sequencing adapter and data was subsequently aligned to fasta files which contained the adapters.
This guarantees the aliugnemnt of coordinates between Illumina, MS and ONT references.
"""  
# Substract length of 24nt long 5' adapter and remove reads which lie in 5' or 3' adapter range (0 <= tRNA_Seq <= 95nts)
# final_native_df = final_native_df.with_columns([
#     (pl.col("chromStart") - 24).alias("chromStart"),
#     (pl.col("chromEnd") - 24).alias("chromEnd"),
# ]).filter((0 <= pl.col("chromStart")) & (pl.col("chromStart") <= 96))

# final_ivt_df = final_ivt_df.with_columns([
#     (pl.col("chromStart") - 24).alias("chromStart"),
#     (pl.col("chromEnd") - 24).alias("chromEnd"),
# ]).filter((0 <= pl.col("chromStart")) & (pl.col("chromStart") <= 96))



# %%
write_bedrmod_polars(final_native_df,list_of_native_dfs_comments[0],output_filename="combined_and_scored.bed")
write_bedrmod_polars(final_ivt_df,list_of_IVT_dfs_comments[0],output_filename="combined_and_scored.bed")
