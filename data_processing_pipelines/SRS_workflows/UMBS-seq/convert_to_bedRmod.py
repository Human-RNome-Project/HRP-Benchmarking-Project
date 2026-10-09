#!/usr/bin/env python3
"""Convert called m5C sites (output of `m5C_calling.py`) into a bedRMod v2 file.

Accepts either a single-sample call table or the long-format intersection
table.  With `--combine-samples` the rows of one site are merged across samples
(depths and unconverted counts are summed), which is what should be used for
the intersection table.

Only the python standard library is required (python >= 3.8).
"""

import argparse
import gzip
import sys

MODS = {"m5C": ("20607", "C"), "Y": ("17802", "U"), "m6A": ("21891", "A")}


def smart_open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def rename_chrom(chrom, style):
    if style == "ucsc":
        if chrom in ("MT", "chrMT"):
            return "chrM"
        return chrom if chrom.startswith("chr") else "chr" + chrom
    if style == "ensembl":
        if chrom in ("chrM", "chrMT", "M"):
            return "MT"
        return chrom[3:] if chrom.startswith("chr") else chrom
    return chrom


def chrom_sort_key(chrom):
    core = chrom[3:] if chrom.startswith("chr") else chrom
    return (0, int(core), "") if core.isdigit() else (1, 0, core)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="called-site table from m5C_calling.py")
    p.add_argument("--output", required=True)
    p.add_argument("--combine-samples", action="store_true",
                   help="merge the rows of one site across samples "
                        "(use for the intersection table)")
    p.add_argument("--mod", default="m5C", choices=sorted(MODS))
    p.add_argument("--name-field", default="short", choices=("short", "id"))
    p.add_argument("--min-coverage", type=int, default=0)
    p.add_argument("--min-frequency", type=float, default=0.0,
                   help="minimum modification percentage, 0-100")
    p.add_argument("--chrom-style", default="asis", choices=("asis", "ucsc", "ensembl"))
    p.add_argument("--organism", default="9606")
    p.add_argument("--assembly", default="GRCh38")
    p.add_argument("--annotation-source", default="Ensembl")
    p.add_argument("--annotation-version", default="110")
    p.add_argument("--platform", default="Illumina NovaSeq X")
    p.add_argument("--workflow",
                   default="UBS-seq (hisat-3n 0.0.3) + binomial m5C calling")
    args = p.parse_args()

    mod_id, mod_base = MODS[args.mod]
    name = mod_id if args.name_field == "id" else args.mod

    sites = {}   # (chrom, pos, strand) -> [depth, unconverted]
    with smart_open(args.input) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for col in ("Chrom", "Pos", "Strand", "Depth", "Unconverted"):
            if col not in idx:
                sys.exit(f"error: column '{col}' not found in {args.input}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (f[idx["Chrom"]], int(f[idx["Pos"]]), f[idx["Strand"]])
            depth, unconv = int(f[idx["Depth"]]), int(f[idx["Unconverted"]])
            if key in sites:
                if not args.combine_samples:
                    continue           # keep the first occurrence
                sites[key][0] += depth
                sites[key][1] += unconv
            else:
                sites[key] = [depth, unconv]

    rows = []
    for (chrom, pos, strand), (depth, unconv) in sites.items():
        frequency = 100.0 * unconv / depth if depth else 0.0
        if depth < args.min_coverage or frequency < args.min_frequency:
            continue
        chrom = rename_chrom(chrom, args.chrom_style)
        start0 = pos - 1
        rows.append((chrom, start0, start0 + 1, name, min(1000, depth), strand,
                     start0, start0 + 1, "0,0,0", depth, f"{frequency:.2f}"))

    rows.sort(key=lambda r: (chrom_sort_key(r[0]), r[1], r[5]))

    with open(args.output, "w") as out:
        out.write("#fileformat=bedRModv2\n")
        out.write(f"#organism={args.organism}\n")
        out.write("#modification_type=RNA\n")
        out.write(f"#modification_names={mod_id}:{args.mod}:{mod_base}\n")
        out.write(f"#assembly={args.assembly}\n")
        out.write(f"#annotation_source={args.annotation_source}\n")
        out.write(f"#annotation_version={args.annotation_version}\n")
        out.write(f"#sequencing_platform={args.platform}\n")
        out.write("#basecalling=\n")
        out.write(f"#bioinformatics_workflow={args.workflow}\n")
        out.write("#experiment=\n")
        out.write("#external_source=\n")
        out.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t"
                  "thickStart\tthickEnd\titemRgb\tcoverage\tfrequency\n")
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    sys.stderr.write(f"{len(rows)} sites written to {args.output}\n")


if __name__ == "__main__":
    main()
