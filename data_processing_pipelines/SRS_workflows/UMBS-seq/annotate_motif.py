#!/usr/bin/env python3
"""Reformat `hisat-3n-table` output and append the trinucleotide context.

Reads the raw hisat-3n-table stream on stdin

    ref  pos  strand  convertedBaseQualities  convertedBaseCount
         unconvertedBaseQualities  unconvertedBaseCount

and writes the per-position table used by `m5C_calling.py`:

    Sample Chrom Pos Strand Converted Unconverted Depth Ratio_conv Ratio_unconv [Motif]

Positions with zero depth are dropped.  Only the python standard library is
required; the reference is read through its `.fai` index, one chromosome at a
time, so memory stays at the size of the largest chromosome.
"""

import argparse
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class Fasta:
    """Minimal random-access FASTA reader backed by the .fai index."""

    def __init__(self, path):
        self.path = path
        self.index = {}
        with open(path + ".fai") as fh:
            for line in fh:
                name, length, offset, line_bases, line_width = line.split("\t")[:5]
                self.index[name] = (int(length), int(offset),
                                    int(line_bases), int(line_width))
        self._name = None
        self._seq = None

    def chrom(self, name):
        if name != self._name:
            if name not in self.index:
                self._name, self._seq = name, None
                return None
            length, offset, line_bases, line_width = self.index[name]
            nlines = (length + line_bases - 1) // line_bases
            with open(self.path, "rb") as fh:
                fh.seek(offset)
                raw = fh.read(nlines * line_width)
            self._seq = raw.replace(b"\n", b"").replace(b"\r", b"")[:length].upper()
            self._name = name
        return self._seq

    def context(self, name, pos, strand):
        """Trinucleotide centred on the 1-based `pos`, on the given strand."""
        seq = self.chrom(name)
        if seq is None or pos < 2 or pos + 1 > len(seq):
            return "NNN"
        tri = seq[pos - 2:pos + 1].decode()
        if strand == "-":
            tri = tri.translate(COMP)[::-1]
        return tri


def classify(tri):
    if len(tri) != 3 or tri[0] != "C":
        return "NA"
    if tri[1] == "G":
        return "CpG"
    if tri[2] == "G":
        return "CHG"
    if "N" in tri:
        return "NA"
    return "CHH"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--ref", help="reference FASTA (must have a .fai index)")
    p.add_argument("--no-motif", action="store_true",
                   help="do not append the Motif column (no --ref needed)")
    p.add_argument("--literal-motif", action="store_true",
                   help="write the trinucleotide itself (e.g. CAG) instead of "
                        "its CpG/CHG/CHH class")
    args = p.parse_args()

    if not args.no_motif and not args.ref:
        sys.exit("error: --ref is required unless --no-motif is given")

    fa = Fasta(args.ref) if not args.no_motif else None

    out = sys.stdout
    cols = ["Sample", "Chrom", "Pos", "Strand", "Converted", "Unconverted",
            "Depth", "Ratio_conv", "Ratio_unconv"]
    if not args.no_motif:
        cols.append("Motif")
    out.write("\t".join(cols) + "\n")

    first = True
    for line in sys.stdin:
        f = line.rstrip("\n").split("\t")
        if first:
            first = False
            if f[0] in ("ref", "#ref"):        # skip the hisat-3n-table header
                continue
        if len(f) < 7:
            continue
        chrom, pos, strand = f[0], int(f[1]), f[2]
        conv, unconv = int(f[4]), int(f[6])
        depth = conv + unconv
        if depth == 0:
            continue
        row = [args.sample, chrom, str(pos), strand, str(conv), str(unconv),
               str(depth), f"{conv / depth:.6g}", f"{unconv / depth:.6g}"]
        if not args.no_motif:
            tri = fa.context(chrom, pos, strand)
            row.append(tri if args.literal_motif else classify(tri))
        out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
