#!/usr/bin/env python3
"""
Normalize a BEDRMod file header to the #[col_list] format.

Converts the BEDRMod v2 tab-separated column-name line (#chrom\t...) to
the Python-list format (#['chrom', ...]) used by this workflow. All other
lines (metadata comments and data rows) are passed through unchanged.
Handles plain and gzip-compressed input; always writes plain output.
"""

import gzip
import sys
import argparse


# Spurious token the raw files inject into the header, both as a trailing
# tab-separated field and as a stray key=value value.
_BOGUS_META = "score_nmod_thresholds"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",  help="Input BED file (plain or .gz)")
    p.add_argument("output", help="Output BED file (plain)")
    p.add_argument(
        "--drop-chroms", action="append", metavar="CHROM", default=[],
        help="Chromosome to exclude from output; repeatable (e.g. --drop-chroms chrY)",
    )
    return p.parse_args()


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def main():
    args = parse_args()
    drop = set(args.drop_chroms)
    converted = False

    with _open(args.input) as fh, open(args.output, "w") as out:
        for line in fh:
            if line.startswith("#"):
                stripped = line.lstrip("#").strip()
                # Tab-format column header (#chrom\t...) → #[list] format.
                # Metadata lines can carry a spurious trailing tab field (an
                # erroneous "score_nmod_thresholds" column); drop it so only the
                # "#key=value" part remains.
                if "\t" in stripped and not stripped.startswith("["):
                    cols = stripped.split("\t")
                    if cols[0] in ("chrom", "chr"):
                        out.write(f"#{cols}\n")        # real column header
                        converted = True
                    else:
                        out.write(f"#{cols[0]}\n")     # metadata: keep key=value
                    continue
                # Blank the same bogus token when it lands as a value, e.g.
                # #experiment=score_nmod_thresholds -> #experiment=
                if "=" in stripped and stripped.split("=", 1)[1] == _BOGUS_META:
                    out.write(f"#{stripped.split('=', 1)[0]}=\n")
                    continue
                out.write(line)
                continue
            if drop and line.split("\t", 1)[0] in drop:
                continue
            out.write(line)

    if not converted:
        print(
            f"[sanitize] note: no tab-format column header found in "
            f"{args.input} (file may already be in #[list] format)",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
