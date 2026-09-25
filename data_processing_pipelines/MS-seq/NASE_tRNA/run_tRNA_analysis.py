#!/usr/bin/env python3
"""Run the Human RNome MS consensus-creation search workflow.

This script reproduces the relevant DecoyDatabase -> NucleicAcidSearchEngine
steps from the TOPPAS workflow. Static search parameters are loaded from the
INI files shipped next to this script. Run-specific inputs and mass tolerances
are supplied on the command line.
"""

from __future__ import annotations

import argparse
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DECOY_INI = SCRIPT_DIR / "decoy_database.ini"
DEFAULT_NASE_INI = SCRIPT_DIR / "NASE.ini"
DEFAULT_CHEBI_MAPPING = SCRIPT_DIR / "ChEBI_ID_RNA_mods_compatible.csv"


def existing_file(value: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"File does not exist: {path}")
    return path


def nonnegative_float(value: str) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Not a number: {value}") from exc
    if number < 0:
        raise argparse.ArgumentTypeError("Mass tolerance must be >= 0")
    return number


def find_executable(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise RuntimeError(
            f"Could not find '{name}' on PATH. Activate/install the OpenMS "
            "environment containing the Human RNome workflow tools first."
        )
    return executable


def run_command(command: list[str]) -> None:
    print(f"+ {shlex.join(command)}", flush=True)
    subprocess.run(command, check=True)


def mzml_stem(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".mzml"):
        return name[:-5]
    return path.stem


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run DecoyDatabase followed by NucleicAcidSearchEngine and write "
            "the resulting bedRMod file. Mass tolerances are in ppm, matching "
            "the Human RNome TOPPAS workflow."
        )
    )
    parser.add_argument("mzml", type=existing_file, help="Input mzML file")
    parser.add_argument("fasta", type=existing_file, help="Input RNA FASTA file")
    parser.add_argument(
        "--precursor-tolerance",
        required=True,
        type=nonnegative_float,
        metavar="PPM",
        help="Precursor mass tolerance in ppm",
    )
    parser.add_argument(
        "--product-tolerance",
        required=True,
        type=nonnegative_float,
        metavar="PPM",
        help="Product/fragment ion mass tolerance in ppm",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for the final .bedrmod file",
    )
    parser.add_argument(
        "--decoy-ini",
        type=existing_file,
        default=DEFAULT_DECOY_INI,
        help=f"DecoyDatabase INI (default: {DEFAULT_DECOY_INI.name})",
    )
    parser.add_argument(
        "--nase-ini",
        type=existing_file,
        default=DEFAULT_NASE_INI,
        help=f"NucleicAcidSearchEngine INI (default: {DEFAULT_NASE_INI.name})",
    )
    parser.add_argument(
        "--chebi-mapping",
        type=existing_file,
        default=DEFAULT_CHEBI_MAPPING,
        help=f"bedRMod ChEBI mapping CSV (default: {DEFAULT_CHEBI_MAPPING.name})",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    # argparse does not apply `type=` to default Path objects, so validate the
    # bundled defaults explicitly as well.
    for label, path in (
        ("DecoyDatabase INI", Path(args.decoy_ini)),
        ("NucleicAcidSearchEngine INI", Path(args.nase_ini)),
        ("ChEBI mapping", Path(args.chebi_mapping)),
    ):
        if not path.is_file():
            raise RuntimeError(f"{label} not found: {path}")

    decoy_database = find_executable("DecoyDatabase")
    nase = find_executable("NucleicAcidSearchEngine")

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    stem = mzml_stem(args.mzml)
    bedrmod_out = output_dir / f"{stem}.bedrmod"

    with tempfile.TemporaryDirectory(prefix="hrp_ms_consensus_") as temp_name:
        temp_dir = Path(temp_name)
        decoy_fasta = temp_dir / "target_decoy.fasta"
        mztab_out = temp_dir / f"{stem}.mzTab"

        run_command(
            [
                decoy_database,
                "-ini",
                str(args.decoy_ini),
                "-in",
                str(args.fasta),
                "-out",
                str(decoy_fasta),
            ]
        )

        run_command(
            [
                nase,
                "-ini",
                str(args.nase_ini),
                "-in",
                str(args.mzml),
                "-database",
                str(decoy_fasta),
                "-out",
                str(mztab_out),
                "-bedrmod_out",
                str(bedrmod_out),
                "-bedrmod_chebi_mapping",
                str(args.chebi_mapping),
                "-precursor:mass_tolerance",
                str(args.precursor_tolerance),
                "-fragment:mass_tolerance",
                str(args.product_tolerance),
            ]
        )

    print(f"bedRMod output: {bedrmod_out}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print(
            f"ERROR: command failed with exit status {exc.returncode}",
            file=sys.stderr,
        )
        raise SystemExit(exc.returncode)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
