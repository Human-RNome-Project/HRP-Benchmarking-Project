from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

from pymol import cmd


GITHUB_ROOT = Path(__file__).resolve().parents[2]
INPUTS_DIR = GITHUB_ROOT / "inputs"
OUTPUTS_DIR = GITHUB_ROOT / "outputs"
OUTDIR = GITHUB_ROOT / "figures" / "panel_h_tiered_rRNA_only_9o3v"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ── Paths ────────────────────────────────────────────────────────────────────
LOCAL_PDB = INPUTS_DIR / "9o3v.cif"
FALLBACK_PDB = Path.home() / "ref" / "9o3v.cif"
PDB_FILE = LOCAL_PDB if LOCAL_PDB.exists() else FALLBACK_PDB

LOCAL_FASTA = INPUTS_DIR / "rcsb_pdb_9O3V.fasta"
FALLBACK_FASTA = Path.home() / "ref" / "rcsb_pdb_9O3V.fasta"
FASTA_FILE = LOCAL_FASTA if LOCAL_FASTA.exists() else FALLBACK_FASTA

TSV_FILE = OUTPUTS_DIR / "tiered_lists" / "tiered_rRNA_only.tsv"

RRNA_18S = "chain S2"
RRNA_28S = "chain L5"
RRNA_5_8S = "chain L8"

MODIFICATION_COLORS = {
    "Am": "#D44F3E",
    "Cm": "#0D3B6E",
    "Gm": "#74B354",
    "Um": "#C47A02",
    "Y": "#F0A202",
    "m5C": "#001427",
    "m6A": "#721817",
}

TRNA_COLOR = "#432818"
NASCENT_COLOR = "#8B1A1A"


def hex_to_rgb(color: str) -> list[float]:
    color = color.lstrip("#")
    return [int(color[i : i + 2], 16) / 255.0 for i in (0, 2, 4)]


def load_structure() -> str:
    objects = cmd.get_names("objects")
    if objects:
        return objects[0]
    cmd.load(str(PDB_FILE), "ribosome")
    return "ribosome"


def parse_chain_sets():
    small = []
    large = []
    for line in FASTA_FILE.read_text().splitlines():
        if not line.startswith(">"):
            continue
        m = re.search(r"\|Chain ([^\[]+)\[auth ([^\]]+)\]\|(.*)$", line)
        if not m:
            continue
        _label_chain, auth_chain, desc = m.groups()
        auth_chain = auth_chain.strip()
        if "40S ribosomal protein" in desc or "Small ribosomal subunit protein" in desc:
            small.append(auth_chain)
        elif (
            "60S ribosomal protein" in desc
            or "Large ribosomal subunit protein" in desc
            or "Ubiquitin-ribosomal protein eL40 fusion protein" in desc
        ):
            large.append(auth_chain)
    return small, large


def sel_from_chains(chains):
    return " or ".join(f"chain {c}" for c in chains)


def read_modifications():
    by_type = defaultdict(list)
    total = 0
    chain_map = {
        "hs_rRNA_18S": "S2",
        "hs_rRNA_28S": "L5",
        "hs_rRNA_5.8S": "L8",
    }
    with TSV_FILE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("in_ref", "").strip().upper() != "TRUE":
                continue
            chrom = row.get("chr", "").strip()
            chain = chain_map.get(chrom)
            if chain is None:
                continue
            mod = (row.get("name") or row.get("ref_mod_type") or "").strip()
            if not mod or mod == "NA":
                mod = (row.get("ref_mod_type") or "").strip()
            if not mod or mod == "NA":
                continue
            resi = int(row["start"]) + 1
            by_type[mod].append((chain, resi))
            total += 1
    return by_type, total


def mod_atom_sel(pairs):
    parts = []
    for chain, resi in pairs:
        parts.append(
            f"(chain {chain} and resi {resi} and not name P+OP1+OP2+OP3+O1P+O2P+O3P)"
        )
    return " or ".join(parts) if parts else "none"


def set_hex_color(name: str, hex_color: str) -> None:
    cmd.set_color(name, hex_to_rgb(hex_color))


def main() -> None:
    obj = load_structure()
    small_proteins, large_proteins = parse_chain_sets()
    mods, total = read_modifications()

    cmd.disable("all")
    cmd.enable(obj)
    cmd.hide("everything", "all")

    cmd.show("cartoon", f"{obj} and polymer")
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("cartoon_side_chain_helper", 0)
    cmd.set("cartoon_sampling", 14)
    cmd.set("cartoon_tube_radius", 0.85)
    cmd.set("cartoon_transparency", 0.80, "polymer")
    cmd.color("gray85", f"{obj} and polymer")

    cmd.color("lightblue", RRNA_18S)
    cmd.set("cartoon_transparency", 0.80, RRNA_18S)
    cmd.color("wheat", RRNA_28S)
    cmd.set("cartoon_transparency", 0.80, RRNA_28S)
    cmd.color("palegreen", RRNA_5_8S)
    cmd.set("cartoon_transparency", 0.80, RRNA_5_8S)

    if small_proteins:
        small_sel = sel_from_chains(small_proteins)
        cmd.color("gray75", small_sel)
        cmd.set("cartoon_transparency", 0.80, small_sel)

    if large_proteins:
        large_sel = sel_from_chains(large_proteins)
        cmd.color("gray80", large_sel)
        cmd.set("cartoon_transparency", 0.80, large_sel)

    for mod, pairs in sorted(mods.items()):
        sel_name = f"mod_{re.sub(r'[^0-9A-Za-z_]+', '_', mod)}"
        atom_sel = mod_atom_sel(pairs)
        if atom_sel == "none":
            continue
        cmd.select(sel_name, atom_sel)
        cmd.show("sticks", sel_name)
        color_name = f"modcol_{re.sub(r'[^0-9A-Za-z_]+', '_', mod)}"
        set_hex_color(color_name, MODIFICATION_COLORS.get(mod, "#FF69B4"))
        cmd.color(color_name, sel_name)
    cmd.set("stick_radius", 0.4, sel_name)
    cmd.set("stick_transparency", 0.0, sel_name)

    set_hex_color("rnome_trna", "#bebcb7")
    set_hex_color("rnome_nascent", "#E2813C")

    cmd.hide("everything", "chain Pt")
    cmd.show("cartoon", "chain Pt")
    cmd.color("rnome_trna", "chain Pt")
    cmd.set("cartoon_transparency", 0.0, "chain Pt")

    cmd.hide("everything", "chain NC")
    cmd.show("spheres", "chain NC")
    cmd.color("rnome_nascent", "chain NC")
    cmd.set("sphere_transparency", 0.0, "chain NC")

    cmd.orient(f"{obj} and polymer")
    cmd.zoom(f"{obj} and polymer", 20)
    cmd.bg_color("white")
    cmd.set("ray_trace_mode", 1)
    cmd.set("antialias", 2)
    
    # Save the PyMOL session image
    output_png = OUTDIR / "panel_h_tiered_rRNA_only_9o3v.png"
    cmd.png(str(output_png), width=2400, height=2400, dpi=300, ray=1)
    
    print(f"tiered_rRNA_only_9o3v prepared ({total} modification sites).")
    print(f"Saved: {output_png}")


main()
