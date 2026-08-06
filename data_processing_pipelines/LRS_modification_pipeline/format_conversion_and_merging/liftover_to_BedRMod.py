# pip install polars numpy argparse tqdm
import polars as pl
import numpy as np
import argparse
from tqdm import tqdm

# Argument Parser
opt_parser = argparse.ArgumentParser()

opt_parser.add_argument(
    "-b",
    "--bed_files",
    dest="bed_files",
    help="Insert an alignment dataframe",
    nargs="+",
)

opt_parser.add_argument(
    "-r",
    "--rna_types",
    dest="rna_types",  # mRNA, tRNA, rRNA
    help="Define the RNA types of the provided bed files",
    nargs="+",
)

opt_parser.add_argument(
    "-t",
    "--technology",
    dest="technology",
    help="Technology in Use -> Illumina, ONT, MS",
)

opt_parser.add_argument(
    "-s",
    "--tRNA_adapter_offset",
    dest="tRNA_adapter_offset",
    help="tRNA adapter offset when using ONT to correct for alignment with adapters (default=-3)",
    default=0,
)

opt_parser.add_argument(
    "-o", "--output_file", dest="output_file", help="Insert an output file to write to"
)


options = opt_parser.parse_args()

bed_files = options.bed_files
rna_types = options.rna_types
tRNA_offset = options.tRNA_adapter_offset
technology = options.technology
output_file = options.output_file

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
    "m2xp7G" : "m2xp7G",

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
    "mA?": "mxA",
    "mxA": "mxA",
    "99990": "mxA",               # placeholder
    "mC?": "mxC",
    "mxC": "mxC",
    "99991": "mxC",               # placeholder
    "mG?": "mxG",
    "mxG": "mxG",
    "99992": "mxG",               # placeholder
    "mU?": "mxU",
    "mxU": "mxU",
    "99993": "mxU",               # placeholder
    "U*": "mxU",

    # --- Non-canonical / experimental placeholders ---
    "xp3Cm": "xp3Cm",
    "xp4U": "xp4U",
    "xp6G": "xp6G",
    "xp7G": "xp7G",
    "m2xp7G": "m2xp7G",
    "f5mA": "f5mA",

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
    "mA?:mxA:A,mC?:mxC:C,mG?:mxG:G,mU?:mxU:U,"
    "mxA:mxA:A,mxC:mxC:C,mxG:mxG:G,mxU:mxU:U,"
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

_MOD_TO_SYMBOL: dict[str, str] = {
    "pmnm5U": "{",
    "m1Am": "œ",
    "pm1Am": "œ",
    "m1Gm": "ε",
    "pm1Gm": "ε",
    "m1Im": "ξ",
    "pm1Im": "ξ",
    "pse2U": "ω",
    "ps2U": "2",
    "m1acp3Y": "α",
    "pm1acp3Y": "α",
    "m1A": "Ѣ",
    "m1G": "K",
    "m1I": "O",
    "pm1I": "O",
    "m1Y": "]",
    "pm1Y": "]",
    "pm1G": "K",
    "pm5Um": "Ħ",
    "pAr(p)": "ʩ",
    "hm5Cm": "¡",
    "phm5Cm": "¡",
    "Am": "ʍ",
    "pAm": "ʍ",
    "Cm": "B",
    "Gm": "＃",
    "Im": "Ш",
    "pIm": "Ш",
    "Ym": "Z",
    "pYm": "Z",
    "Um": "J",
    "mcmo5Um": "Ͽ",
    "pmcmo5Um": "Ͽ",
    "Ar(p)": "ʩ",
    "Gr(p)": "ℑ",
    "pGr(p)": "ℑ",
    "pGr": "Ѥ",
    "N2'3'cp": "Ғ",
    "m2,8A": "±",
    "pm2,8A": "±",
    "msms2i6A": "£",
    "pmsms2i6A": "£",
    "ges2U": "Γ",
    "pges2U": "Γ",
    "k2C": "}",
    "pk2C": "}",
    "m2A": "ɿ",
    "pm2A": "ɿ",
    "ms2ct6A": "ÿ",
    "pms2ct6A": "ÿ",
    "ms2io6A": "≠",
    "pms2io6A": "≠",
    "ms2hn6A": "≈",
    "pms2hn6A": "≈",
    "pms2i6A": "*",
    "ms2i6A": "*",
    "ms2m6A": "∞",
    "pms2m6A": "∞",
    "ms2t6A": "[",
    "pms2t6A": "[",
    "se2U": "ω",
    "s2Um": "∏",
    "ps2Um": "∏",
    "s2C": "ʤ",
    "ps2C": "ʤ",
    "s2U": "2",
    "pm2G": "L",
    "m3Um": "σ",
    "pm3Um": "σ",
    "acp3D": "Ð",
    "pacp3D": "Ð",
    "acp3Um": "‡",
    "acp3Y": "Þ",
    "pacp3Y": "Þ",
    "acp3U": "X",
    "pacp3U": "X",
    "pm3C": "Щ",
    "m3C": "Щ",
    "m3Y": "Ƒ",
    "pm3Y": "Ƒ",
    "m3U": "δ",
    "pm3U": "δ",
    "imG-14": "†",
    "pimG-14": "†",
    "pm4C": "ν",
    "s4U": "4",
    "ps4U": "4",
    "pm4Cm": "λ",
    "CoApN": "♠",
    "acCoApN": "♣",
    "malonyl-CoApN": "♥",
    "succinyl-CoApN": "♦",
    "ppN": "ϒ",
    "5'-OH-N": "ɀ",
    "NADpN": "Ξ",
    "pppN": "ϖ",
    "m5Cm": "τ",
    "pm5Cm": "τ",
    "m5Um": "Ħ",
    "pD": "D",
    "pmcm5s2U": "3",
    "mchm5Um": "b",
    "pmchm5Um": "b",
    "mchm5U": "ɮ",
    "pmchm5U": "ɮ",
    "pcmo5U": "V",
    "phm5C": "Ƣ",
    "inm5Um": "Ю",
    "pinm5Um": "Ю",
    "inm5s2U": "Ɲ",
    "pinm5s2U": "Ɲ",
    "inm5U": "¾",
    "pinm5U": "¾",
    "nm5ges2U": "Δ",
    "pnm5ges2U": "Δ",
    "nm5se2U": "π",
    "pnm5se2U": "π",
    "nm5s2U": "∫",
    "pnm5s2U": "∫",
    "nm5U": "∪",
    "pnm5U": "∪",
    "nchm5U": "r",
    "pnchm5U": "r",
    "ncm5Um": "~",
    "pncm5Um": "~",
    "ncm5s2U": "l",
    "pncm5s2U": "l",
    "ncm5U": "&",
    "pncm5U": "&",
    "chm5U": "≥",
    "pchm5U": "≥",
    "cm5s2U": "℘",
    "pcm5s2U": "℘",
    "cmnm5Um": ")",
    "pcmnm5Um": ")",
    "cmnm5ges2U": "f",
    "pcmnm5ges2U": "f",
    "cmnm5se2U": "⊥",
    "pcmnm5se2U": "⊥",
    "cmnm5s2U": "$",
    "pcmnm5s2U": "$",
    "cmnm5U": "!",
    "pcmnm5U": "!",
    "cm5U": "◊",
    "pcm5U": "◊",
    "cnm5U": "Ѷ",
    "pcnm5U": "Ѷ",
    "f5Cm": "°",
    "pf5Cm": "°",
    "f5se2U": None,
    "f5s2U": None,
    "f5Um": None,
    "f5C": ">",
    "pf5C": ">",
    "f5U": None,
    "ho5C": "Ç",
    "pho5C": "Ç",
    "hm5C": "Ƣ",
    "ho5U": "∝",
    "pho5U": "∝",
    "mcm5Um": "∩",
    "pmcm5Um": "∩",
    "mcm5s2U": "3",
    "mcm5U": "1",
    "pmcm5U": "1",
    "mo5U": "5",
    "pmo5U": "5",
    "m5s2U": "F",
    "pm5s2U": "F",
    "mnm5ges2U": "h",
    "pmnm5ges2U": "h",
    "mnm5se2U": "≅",
    "pmnm5se2U": "≅",
    "mnm5s2U": "S",
    "pmnm5s2U": "S",
    "mnm5U": "{",
    "m5C": "?",
    "pm5C": "?",
    "m5D": "ρ",
    "pm5D": "ρ",
    "m5U": "T",
    "pm5U": "T",
    "tm5s2U": "ƕ",
    "ptm5s2U": "ƕ",
    "tm5U": "ʭ",
    "ptm5U": "ʭ",
    "pm6,6A": "ζ",
    "yW-86": "¥",
    "pyW-86": "¥",
    "yW-72": "Ω",
    "yW-58": "⇑",
    "pyW-58": "⇑",
    "pyW-72": "Ω",
    "preQ1base": "∇",
    "preQ1": "∉",
    "ppreQ1": "∉",
    "preQ0base": "ψ",
    "preQ0": "φ",
    "ppreQ0": "φ",
    "m7G": "7",
    "pm7G": "7",
    "m8A": "â",
    "pm8A": "â",
    "pm1A": "Ѣ",
    "pac4C": "M",
    "m2Gm": "γ",
    "pm2Gm": "γ",
    "m2,7Gm": "æ",
    "pm2,7Gm": "∨",
    "m2,7G": "∨",
    "m2,7GpppN": "®",
    "pm2,7G": "æ",
    "m2,2Gm": "|",
    "pm2,2Gm": "|",
    "m2,2,7G": "∠",
    "m2,2,7GpppN": "¶",
    "pm2,2,7G": "∠",
    "m2,2G": "R",
    "pm2,2G": "R",
    "m2G": "L",
    "m4Cm": "λ",
    "m4,4Cm": "β",
    "pm4,4Cm": "β",
    "m4,4C": "μ",
    "ac4Cm": "ℵ",
    "ac4C": "M",
    "m4C": "ν",
    "m6Am": "χ",
    "pm6Am": "χ",
    "m6,6Am": "η",
    "pm6,6Am": "η",
    "m6,6A": "ζ",
    "io6A": "Ỽ",
    "pio6A": "Ỽ",
    "ac6A": "⇓",
    "pac6A": "⇓",
    "f6A": "Ϩ",
    "pf6A": "Ϩ",
    "g6A": "≡",
    "pg6A": "≡",
    "hm6A": "Ϫ",
    "phm6A": "Ϫ",
    "hn6A": "√",
    "phn6A": "√",
    "pi6A": "Ч",
    "i6A": "Ч",
    "m6t6A": "E",
    "pm6t6A": "E",
    "m6ApppppN": "Ϛ",
    "m6AppppN": "Ϙ",
    "m6ApppN": "Ϯ",
    "m6A": "Ж",
    "pm6A": "Ж",
    "t6A": "6",
    "pt6A": "6",
    "m7GppppN": "Ϟ",
    "m7GpppN": "©",
    "pCm": "B",
    "pGm": "＃",
    "pUm": "J",
    "MU3": "亱",
    "pac4Cm": "ℵ",
    "pm4,4C": "μ",
    "A": "A",
    "ApppppN": "Ϭ",
    "AppppN": "Ĳ",
    "ApppN": "Ϩ",
    "pAp": "Ợ",
    "ppA": None,
    "pA": "A",
    "pA2'3'cp": "Ҕ",
    "pppA": "ƿ",
    "C+": "¿",
    "pC+": "¿",
    "mmpN": "Ѽ",
    "mpN": "Ѿ",
    "Nm": "Ѭ",
    "N": ".",
    "G+": "(",
    "pG+": "(",
    "m7GpppNm": None,
    "m7GpppNmNm": None,
    "m7GpppNNm": None,
    "GpppNmN": None,
    "GpppNm": None,
    "ct6A": "e",
    "pct6A": "e",
    "C": "C",
    "pC": "C",
    "pC2'3'cp": "Ҏ",
    "D": "D",
    "oQ": "ς",
    "poQ": "ς",
    "galQ": "9",
    "pgalQ": "9",
    "mpppN": "§",
    "gluQ": "⊄",
    "pgluQ": "⊄",
    "pG2'3'cp": "Ґ",
    "G": "G",
    "pG(pN)": "Ʉ",
    "GpppN": "ϑ",
    "pGp": "Ⱥ",
    "ppG": "Ƈ",
    "pG": "G",
    "pppG": "Ȝ",
    "ht6A": "«",
    "pht6A": "«",
    "OHyW": "⊆",
    "pOHyW": "⊆",
    "I": "I",
    "pI": "I",
    "imG2": "⊇",
    "pimG2": "⊇",
    "manQ": "8",
    "pmanQ": "8",
    "OHyWy": "y",
    "pOHyWy": "y",
    "mimG": "∑",
    "pmimG": "∑",
    "o2yW": "W",
    "po2yW": "W",
    "Y": "P",
    "pY": "P",
    "Qbase": "∴",
    "Q": "Q",
    "pQ": "Q",
    "OHyWx": "š",
    "pOHyWx": "š",
    "pN": "m",
    "xX": "Ѯ",
    "xA": "H",
    "xC": "Ѵ",
    "xG": "ʆ",
    "xU": "N",
    "Xm": "Ѩ",
    "U": "U",
    "cmo5U": "V",
    "mcmo5U": "υ",
    "pmcmo5U": "υ",
    "pU": "U",
    "pU2'3'cp": "Ҋ",
    "yW": "Y",
    "pyW": "Y",
    "pyyW": "ϓ",
    "imG": "€",
    "pimG": "€",
    "mxA": "mxA",
    "mxC": "mxC",
    "mxG": "mxG",
    "mxU": "mxU",
    "U*": "U*",
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
    "T": "T",
    "m7A": "m7A",
    "m6G": "m6G",
    "ac7G": "ac7G",
    "xp3Cm": "xp3Cm",
    "xp4U": "xp4U",
    "xp6G": "xp6G",
    "xp7G": "xp7G",
    "m2xp7G": "m2xp7G",
    "f5mA" : "f5mA",
}

# Dictionary 2: modification name → numeric code
_MOD_TO_ID: dict[str, str] = {
    "pmnm5U": "2000511551U",
    "m1Am": "2000000901A",
    "pm1Am": "2000901551A",
    "m1Gm": "2000000901G",
    "pm1Gm": "2000901551G",
    "m1Im": "2000009019A",
    "pm1Im": "2009019551A",
    "pse2U": "2000020551U",
    "ps2U": "2000002551U",
    "m1acp3Y": "2000001309U",
    "pm1acp3Y": "2001309551U",
    "m1A": "2000000001A",
    "m1G": "2000000001G",
    "m1I": "2000000019A",
    "pm1I": "2000019551A",
    "m1Y": "2000000019U",
    "pm1Y": "2000019551U",
    "pm1G": "2000001551G",
    "pm5Um": "2000905551U",
    "pAr(p)": "2900255551A",
    "hm5Cm": "2000009051C",
    "phm5Cm": "2009051551C",
    "Am": "2000000090A",
    "pAm": "2000090551A",
    "Cm": "2000000090C",
    "Gm": "2000000090G",
    "Im": "2000000909A",
    "pIm": "2000909551A",
    "Ym": "2000000909U",
    "pYm": "2000909551U",
    "Um": "2000000090U",
    "mcmo5Um": "2000090503U",
    "pmcmo5Um": "2090503551U",
    "Ar(p)": "2000000000A",
    "Gr(p)": "2000000900G",
    "pGr(p)": "2900255551G",
    "pGr": "2000900551G",
    "N2'3'cp": "2000003377N",
    "m2,8A": "2000000028A",
    "pm2,8A": "2000028551A",
    "msms2i6A": "2000021161A",
    "pmsms2i6A": "2021161551A",
    "ges2U": "2000000021U",
    "pges2U": "2000021551U",
    "k2C": "2000000021C",
    "pk2C": "2000021551C",
    "m2A": "2000000002A",
    "pm2A": "2000002551A",
    "ms2ct6A": "2000002164A",
    "pms2ct6A": "2002164551A",
    "ms2io6A": "2000002160A",
    "pms2io6A": "2002160551A",
    "ms2hn6A": "2000002163A",
    "pms2hn6A": "2002163551A",
    "pms2i6A": "2002161551A",
    "ms2i6A": "2000002161A",
    "ms2m6A": "2000000621A",
    "pms2m6A": "2000621551A",
    "ms2t6A": "2000002162A",
    "pms2t6A": "2002162551A",
    "se2U": "2000000020U",
    "s2Um": "2000000902U",
    "ps2Um": "2000902551U",
    "s2C": "2000000002C",
    "ps2C": "2000002551C",
    "s2U": "2000000002U",
    "pm2G": "2000002551G",
    "m3Um": "2000000903U",
    "pm3Um": "2000903551U",
    "acp3D": "2000000308U",
    "pacp3D": "2000308551U",
    "acp3Um": None,
    "acp3Y": "2000000309U",
    "pacp3Y": "2000309551U",
    "acp3U": "2000000030U",
    "pacp3U": "2000030551U",
    "pm3C": "2000003551C",
    "m3C": "2000000003C",
    "m3Y": "2000000039U",
    "pm3Y": "2000039551U",
    "m3U": "2000000003U",
    "pm3U": "2000003551U",
    "imG-14": "2000000004G",
    "pimG-14": "2000004551G",
    "pm4C": "2000004551C",
    "s4U": "2000000074U",
    "ps4U": "2000074551U",
    "pm4Cm": "2000904551C",
    "CoApN": "2000000455N",
    "acCoApN": "2000004155N",
    "malonyl-CoApN": "2000004255N",
    "succinyl-CoApN": "2000004355N",
    "ppN": "2000000552N",
    "5'-OH-N": "2000000550N",
    "NADpN": "2000000255N",
    "pppN": "2000000553N",
    "m5Cm": "2000000905C",
    "pm5Cm": "2000905551C",
    "m5Um": "2000000905U",
    "pD": "2000008551U",
    "pmcm5s2U": "2002521551U",
    "mchm5Um": "2000090522U",
    "pmchm5Um": "2090522551U",
    "mchm5U": "2000000522U",
    "pmchm5U": "2000522551U",
    "pcmo5U": "2000502551U",
    "phm5C": "2000051551C",
    "inm5Um": "2000090583U",
    "pinm5Um": "2090583551U",
    "inm5s2U": "2000002583U",
    "pinm5s2U": "2002583551U",
    "inm5U": "2000000583U",
    "pinm5U": "2000583551U",
    "nm5ges2U": "2000021510U",
    "pnm5ges2U": "2021510551U",
    "nm5se2U": "2000020510U",
    "pnm5se2U": "2020510551U",
    "nm5s2U": "2000002510U",
    "pnm5s2U": "2002510551U",
    "nm5U": "2000000510C",
    "pnm5U": "2000510551U",
    "nchm5U": "2000000531U",
    "pnchm5U": "2000531551U",
    "ncm5Um": "2000009053U",
    "pncm5Um": "2009053551U",
    "ncm5s2U": "2000000253U",
    "pncm5s2U": "2000253551U",
    "ncm5U": "2000000053U",
    "pncm5U": "2000053551U",
    "chm5U": "2000000520U",
    "pchm5U": "2000520551U",
    "cm5s2U": "2000002540U",
    "pcm5s2U": "2002540551U",
    "cmnm5Um": "2000009051U",
    "pcmnm5Um": "2009051551U",
    "cmnm5ges2U": "2000002151U",
    "pcmnm5ges2U": "2002151551U",
    "cmnm5se2U": "2000002051U",
    "pcmnm5se2U": "2002051551U",
    "cmnm5s2U": "2000000251U",
    "pcmnm5s2U": "2000251551U",
    "cmnm5U": "2000000051U",
    "pcmnm5U": "2000051551U",
    "cm5U": "2000000052U",
    "pcm5U": "2000052551U",
    "cnm5U": "2000000055U",
    "pcnm5U": "2000055551U",
    "f5Cm": "2000009071C",
    "pf5Cm": "2009071551C",
    "f5se2U": None,
    "f5s2U": None,
    "f5Um": None,
    "f5C": "2000000071C",
    "pf5C": "2000071551C",
    "f5U": None,
    "ho5C": "2000000050C",
    "pho5C": "2000050551C",
    "hm5C": "2000000051C",
    "ho5U": "2000000050U",
    "pho5U": "2000050551U",
    "mcm5Um": "2000090521U",
    "pmcm5Um": "2090521551U",
    "mcm5s2U": "2000002521U",
    "mcm5U": "2000000521U",
    "pmcm5U": "2000521551U",
    "mo5U": "2000000501U",
    "pmo5U": "2000501551U",
    "m5s2U": "2000000025U",
    "pm5s2U": "2000025551U",
    "mnm5ges2U": "2000021511U",
    "pmnm5ges2U": "2021511551U",
    "mnm5se2U": "2000020511U",
    "pmnm5se2U": "2020511551U",
    "mnm5s2U": "2000002511U",
    "pmnm5s2U": "2002511551U",
    "mnm5U": "2000000511U",
    "m5C": "2000000005C",
    "pm5C": "2000005551C",
    "m5D": "2000000058U",
    "pm5D": "2000058551U",
    "m5U": "2000000005U",
    "pm5U": "2000005551U",
    "tm5s2U": "2000000254U",
    "ptm5s2U": "2000254551U",
    "tm5U": "2000000054U",
    "ptm5U": "2000054551U",
    "pm6,6A": "2000066551A",
    "yW-86": "2000000047C",
    "pyW-86": "2000047551G",
    "yW-72": "2000000347G",
    "yW-58": "2000000348G",
    "pyW-58": "2000348551G",
    "pyW-72": "2000347551G",
    "preQ1base": "2000101000G",
    "preQ1": "2000000101G",
    "ppreQ1": "2000101551G",
    "preQ0base": "200100000G",
    "preQ0": "2000000100G",
    "ppreQ0": "2000100551G",
    "m7G": "2000000007G",
    "pm7G": "2000007551G",
    "m8A": "2000000008A",
    "pm8A": "2000008551A",
    "pm1A": "2000001551A",
    "pac4C": "2000042551C",
    "m2Gm": "2000000902G",
    "pm2Gm": "2000902551G",
    "m2,7Gm": "2000009027G",
    "pm2,7Gm": "2009027551G",
    "m2,7G": "2000000027G",
    "m2,7GpppN": "2000279553N",
    "pm2,7G": "2000004553N",
    "m2,2Gm": "2000009022G",
    "pm2,2Gm": "2009022551G",
    "m2,2,7G": "2000000227G",
    "m2,2,7GpppN": "2002279553N",
    "pm2,2,7G": "2000227551G",
    "m2,2G": "2000000022G",
    "pm2,2G": "2000022551G",
    "m2G": "2000000002G",
    "m4Cm": "2000000904C",
    "m4,4Cm": "2000009044C",
    "pm4,4Cm": "2009044551C",
    "m4,4C": "2000000044C",
    "ac4Cm": "2000009042C",
    "ac4C": "2000000042C",
    "m4C": "2000000004C",
    "m6Am": "2000000906A",
    "pm6Am": "2000906551A",
    "m6,6Am": "2000009066A",
    "pm6,6Am": "2009066551A",
    "m6,6A": "2000000066A",
    "io6A": "2000000060A",
    "pio6A": "2000060551A",
    "ac6A": "2000000064A",
    "pac6A": "2000064551A",
    "f6A": "2000000067A",
    "pf6A": "2000067551A",
    "g6A": "2000000065A",
    "pg6A": "2000065551A",
    "hm6A": "2000000068A",
    "phm6A": "2000068551A",
    "hn6A": "2000000063A",
    "phn6A": "2000063551A",
    "pi6A": "2000061551A",
    "i6A": "2000000061A",
    "m6t6A": "2000000662A",
    "pm6t6A": "2000662551A",
    "m6ApppppN": "2000064555N",
    "m6AppppN": "2000064554N",
    "m6ApppN": "2000064553N",
    "m6A": "2000000006A",
    "pm6A": "2000006551A",
    "t6A": "2000000062A",
    "pt6A": "2000062551A",
    "m7GppppN": "2000079554N",
    "m7GpppN": "2000079553N",
    "pCm": "2000090551C",
    "pGm": "2000090551G",
    "pUm": "2000090551U",
    "MU3": "40000035551U",
    "pac4Cm": "2009042551C",
    "pm4,4C": "2000044551C",
    "A": "A",
    "ApppppN": "2000004555N",
    "AppppN": "2000004554N",
    "ApppN": "2000004553N",
    "pAp": "2000033551A",
    "ppA": "2000000552A",
    "pA": "A",
    "pA2'3'cp": "2000003377A",
    "pppA": "2000000553A",
    "C+": "2000000020C",
    "pC+": "2000020551C",
    "mmpN": "2000002551N",
    "mpN": "2000001551N",
    "Nm": "2000000090N",
    "N": "2000000000N",
    "G+": "2000000103G",
    "pG+": "2000103551G",
    "m7GpppNm": None,
    "m7GpppNmNm": None,
    "m7GpppNNm": None,
    "GpppNmN": None,
    "GpppNm": None,
    "ct6A": "2000000069A",
    "pct6A": "2000069551A",
    "C": "C",
    "pC": "C",
    "pC2'3'cp": "2000003377C",
    "D": "2000000008U",
    "oQ": "2000000102G",
    "poQ": "2000102551G",
    "galQ": "2000000104G",
    "pgalQ": "2000104551G",
    "mpppN": "2000001553N",
    "gluQ": "2000000105G",
    "pgluQ": "2000105551G",
    "pG2'3'cp": "2000003377G",
    "G": "G",
    "pG(pN)": "2000000551G2000000511N",
    "GpppN": "2000009553N",
    "pGp": "2000033551G",
    "ppG": "2000000552G",
    "pG": "G",
    "pppG": "2000000553G",
    "ht6A": "2000002165A",
    "pht6A": "2002165551A",
    "OHyW": "2000034830G",
    "pOHyW": "2034830551G",
    "I": "2000000009A",
    "pI": "2000009551A",
    "imG2": "2000000042G",
    "pimG2": "2000042551G",
    "manQ": "2000000106G",
    "pmanQ": "2000106551G",
    "OHyWy": "2000003480G",
    "pOHyWy": "2003480551G",
    "mimG": "2000000342G",
    "pmimG": "2000342551G",
    "o2yW": "2000034832G",
    "po2yW": "2034832551G",
    "Y": "2000000009U",
    "pY": "2000009551U",
    "Qbase": "2000010000G",
    "Q": "2000000010G",
    "pQ": "2000042551G",
    "OHyWx": "2000003470G",
    "pOHyWx": "2003470551G",
    "pN": "2000000551N",
    "xX": "2000000000X",
    "xA": "2000000999A",
    "xC": "2000000999C",
    "xG": "2000000999G",
    "xU": "2000000999U",
    "Xm": "2000000090X",
    "U": "U",
    "cmo5U": "2000000502U",
    "mcmo5U": "2000000503U",
    "pmcmo5U": "2000503551U",
    "pU": "U",
    "pU2'3'cp": "2003377551U",
    "yW": "2000003483G",
    "pyW": "2003483551G",
    "pyyW": "2034831551G",
    "imG": "2000000034G",
    "pimG": "2000034551G",
    "mxA": "mxA",
    "mxC": "mxC",
    "mxG": "mxG",
    "mxU": "mxU",
    "U*": "U*",
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
    "T": "T",
    "m7A": "m7A",
    "m6G": "m6G",
    "ac7G": "ac7G",
    "xp3Cm": "xp3Cm",
    "xp4U": "xp4U",
    "xp6G": "xp6G",
    "xp7G": "xp7G",
    "m2xp7G": "m2xp7G",
    "f5mA" : "f5mA",
}


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

    df = df.select(df.columns[:12])

    n_cols = len(df.columns)
    if n_cols in COLUMN_MAP:
        df = df.rename(dict(zip(df.columns, COLUMN_MAP[n_cols])))
    else:
        print(f"Warning: unexpected column count ({n_cols}); columns left unnamed.")

    if df["frequency"].max() <= 1.0:
        df = df.with_columns((pl.col("frequency") * 100).alias("frequency"))
    print(df["name"])
    print(df["score"])
    df = df.select(
        [
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
        ]
    ).cast({
    "chrom":       pl.Utf8,
    "chromStart":  pl.Int64,
    "chromEnd":    pl.Int64,
    "name":        pl.Utf8,
    "score":       pl.Float64,      # BED spec: 0–1000
    "strand":      pl.Categorical, # only "+", "-", "."
    "thickStart":  pl.Int64,
    "thickEnd":    pl.Int64,
    "itemRgb":     pl.Utf8,       # e.g. "255,0,0"
    "coverage":    pl.Float64,
    "frequency":   pl.Float64,
    })



    if modtype:
        df = df.with_columns(pl.lit(modtype).alias("name"))

    print(df.head())
    df = align_modification_names_polars(df)
    df = align_chromosome_names_polars(df)
    print("Chromosomes and modifications aligned")

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


def main() -> None:
    """
    Entry point for processing and merging bedRMod files.

    Iterates over all input bed files paired with their RNA type labels, loads each
    file using load_bedrmod_polars, optionally applies a coordinate offset for ONT tRNA
    files, then concatenates all DataFrames vertically. Appends single-letter modification
    codes and modification IDs as new columns, merges the per-file comment headers into a
    single combined header, and writes the final DataFrame to the specified output file.

    Reads from module-level variables:
        bed_files (list[str]): Paths to input bedRMod files.
        rna_types (list[str]): RNA type label for each corresponding input file.
        technology (str): Sequencing technology string (e.g. 'ONT', 'Illumina').
        tRNA_offset (int): Coordinate offset applied to tRNA files when technology is 'ONT'.
        output_file (str): Path for the merged output bedRMod file.

    Returns:
        None
    """
    list_of_processed_dataframes = []
    list_of_processed_comments = []
    for bed_file, rna_type in tqdm(
        zip(bed_files, rna_types), total=len(bed_files), desc="Processing bedrmod files"
    ):
        bed_df, bed_comments = load_bedrmod_polars(bed_file)
        if rna_type == "tRNA" and technology == "ONT":
            bed_df = shift_coordinates(bed_df, offset=tRNA_offset)
        list_of_processed_dataframes.append(bed_df)
        list_of_processed_comments.append(bed_comments)
    processed_bed_df = pl.concat(list_of_processed_dataframes, how="vertical")
    for i in processed_bed_df["name"]:
        try:
            x = _MOD_TO_SYMBOL[i]
            y = _MOD_TO_ID[i]
        except KeyError:
            print(i)
    processed_bed_df = processed_bed_df.with_columns(
        [pl.col("name").replace_strict(_MOD_TO_SYMBOL).alias("single_letter_code")]
    )
    processed_bed_df = processed_bed_df.with_columns(
        pl.col("name").replace_strict(_MOD_TO_ID).alias("mod_id")
    )
    processed_comments = concatenate_bed_comments(
        list_of_processed_comments, rna_types, technology
    )
    write_bedrmod_polars(processed_bed_df, processed_comments, output_file)


if __name__ == "__main__":
    main()
