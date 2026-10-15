#!/usr/bin/env python3
"""Convert a pseudoU-BIDseq `filter_sites` table into a bedRMod v2 file.

The input is `workspace/filter_sites/genome.tsv.gz` (genomic coordinates) or
`workspace/filter_sites/genes.tsv.gz` (coordinates local to a sequence of the
`genes` reference, e.g. the 45S pre-rRNA unit).  Local coordinates are lifted
onto the assembly with `--lift`.

Only the python standard library is required (python >= 3.8).
"""

import argparse
import gzip
import sys

# short name -> (identifier used in the #modification_names header, modified base)
MODS = {
    "Y": ("17802", "U"),
    "m5C": ("20607", "C"),
    "m6A": ("21891", "A"),
}


def smart_open(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_lift(values):
    """'NAME=CHROM:START:STRAND' -> {NAME: (CHROM, START, STRAND)}"""
    lift = {}
    for item in values or []:
        try:
            name, spec = item.split("=", 1)
            chrom, start, strand = spec.split(":")
            lift[name] = (chrom, int(start), strand)
        except ValueError:
            sys.exit(f"error: cannot parse --lift '{item}' "
                     "(expected NAME=CHROM:START:STRAND)")
        if strand not in ("+", "-"):
            sys.exit(f"error: strand of --lift '{item}' must be + or -")
        if strand == "-":
            sys.exit("error: --lift currently supports only sequences extracted "
                     "from the + strand of the assembly")
    return lift


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
    if core.isdigit():
        return (0, int(core), "")
    return (1, 0, core)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="filter_sites/{genes,genome}.tsv.gz")
    p.add_argument("--output", required=True, help="output bedRMod file")
    p.add_argument("--treated-samples", required=True,
                   help="comma-separated sample names of the treated libraries; "
                        "their depths are summed into the coverage column")
    p.add_argument("--group", default=None,
                   help="group name used in the configuration file "
                        "(default: inferred from the *_passed column)")
    p.add_argument("--mod", default="Y", choices=sorted(MODS),
                   help="modification recorded in the file (default: Y)")
    p.add_argument("--name-field", default="short", choices=("short", "id"),
                   help="write the short name ('Y') or the identifier ('17802') "
                        "into the bedRMod name column (default: short)")
    p.add_argument("--passed-only", action="store_true",
                   help="keep only sites with {group}_passed == 1")
    p.add_argument("--min-coverage", type=int, default=0)
    p.add_argument("--min-frequency", type=float, default=0.0,
                   help="minimum modification percentage, 0-100")
    p.add_argument("--lift", action="append", metavar="NAME=CHROM:START:STRAND",
                   help="lift local coordinates of reference sequence NAME onto "
                        "the assembly; may be given several times")
    p.add_argument("--chrom-style", default="asis",
                   choices=("asis", "ucsc", "ensembl"))
    p.add_argument("--organism", default="9606")
    p.add_argument("--assembly", default="GRCh38")
    p.add_argument("--annotation-source", default="Ensembl")
    p.add_argument("--annotation-version", default="110")
    p.add_argument("--platform", default="Illumina NovaSeq X")
    p.add_argument("--workflow", default="pseudoU-BIDseq v2.0 (docker://y9ch/bidseq)")
    args = p.parse_args()

    lift = parse_lift(args.lift)
    treated = [s.strip() for s in args.treated_samples.split(",") if s.strip()]
    mod_id, mod_base = MODS[args.mod]
    name = mod_id if args.name_field == "id" else args.mod

    rows = []
    with smart_open(args.input) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}

        group = args.group
        if group is None:
            passed = [c[:-len("_passed")] for c in header if c.endswith("_passed")]
            if len(passed) != 1:
                sys.exit("error: cannot infer --group from the header "
                         f"(candidates: {passed or 'none'})")
            group = passed[0]

        for col in [f"{group}_fraction", f"{group}_passed"] + \
                   [f"{s}_depth" for s in treated]:
            if col not in idx:
                sys.exit(f"error: column '{col}' not found in {args.input}\n"
                         f"       available: {', '.join(header)}")

        for line in fh:
            f = line.rstrip("\n").split("\t")
            if args.passed_only and f[idx[f"{group}_passed"]] not in ("1", "1.0", "True"):
                continue

            chrom, pos, strand = f[0], int(f[1]), f[2]
            if chrom in lift:
                target, start, _ = lift[chrom]
                chrom, pos = target, start + pos - 1

            chrom = rename_chrom(chrom, args.chrom_style)
            coverage = sum(int(f[idx[f"{s}_depth"]]) for s in treated)
            frequency = float(f[idx[f"{group}_fraction"]]) * 100.0
            if coverage < args.min_coverage or frequency < args.min_frequency:
                continue

            start0 = pos - 1
            rows.append((chrom, start0, start0 + 1, name, min(1000, coverage),
                         strand, start0, start0 + 1, "0,0,0", coverage,
                         f"{frequency:.2f}"))

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
