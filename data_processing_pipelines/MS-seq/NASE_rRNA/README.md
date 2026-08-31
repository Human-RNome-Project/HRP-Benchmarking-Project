# NASE rRNA
This folder contains the rRNA-specific wrapper for the Human RNome MS-seq search workflow. It reproduces the relevant OpenMS DecoyDatabase to NucleicAcidSearchEngine
steps on an `.mzML` file and an RNA FASTA reference, then writes a fragments containing `bedRMod` output file suitable for downstream consensus creation.

The script in this folder is the rRNA analogue of the tRNA workflow and is designed to work with the bundled configuration files in this directory.


## Installation
A conda package will be available soon. For now use the provided [`Dockerfile`](Dockerfile). 

**Build the docker container:**
```bash
    docker build -t hrp-ms-seq:latest .
```

## Input data
The workflow expects:

- an input `.mzML` file containing your MS-seq data
- a FASTA sequence file for the target RNA species 
- configuration files in this folder:
  - `decoy_database.ini`
  - `NASE.ini`
  - `ChEBI_ID_RNA_mods_compatible.csv`

The default rRNA reference in this folder is:

- `hs_rRNAs.fa`

### Usage

**main pipeline**
```bash
docker run --rm -v "$PWD":/data hrp-ms-seq:latest rrna \
    <input.mzML> \
    <input.fasta> \
    --precursor-tolerance 10 \
    --product-tolerance 20 \
    --output-dir <./results> \
    --decoy-ini /data/decoy_database.ini \
    --nase-ini /data/NASE.ini \
    --chebi-mapping /data/ChEBI_ID_RNA_mods_compatible.csv
```

### Arguments

| Argument | Required | Description |
|---|---|---|
| `mzml` | Yes | Input mzML file containing the MS data. |
| `fasta` | Yes | RNA FASTA reference used for the search. |
| `--precursor-tolerance` | Yes | Precursor mass tolerance in ppm. |
| `--product-tolerance` | Yes | Product/fragment ion mass tolerance in ppm. |
| `--output-dir` | Yes | Directory where the result `.bedrmod` file will be written. |
| `--decoy-ini` | No | Custom DecoyDatabase INI file. Defaults to `decoy_database.ini`. |
| `--nase-ini` | No | Custom NucleicAcidSearchEngine INI file. Defaults to `NASE.ini`. |
| `--chebi-mapping` | No | ChEBI mapping CSV for `bedRMod` output. Defaults to `ChEBI_ID_RNA_mods_compatible.csv`. |

## Output
The script writes a fragment containing `bedRMod` file to the requested output directory.

The output file name is based on the input mzML stem:

```text
<output_dir>/<mzml_stem>.bedrmod
```

This output file is intended to be used as input for the downstream consensus creation workflow in the parent [MS-seq README](../README.md).


## Typical workflow
1. Prepare the mzML file.
2. Select the appropriate rRNA FASTA reference.
3. Run the wrapper with your desired precursor/product tolerances.
4. Use the resulting `.bedrmod` file as input for the downstream filtering/consensus pipeline in [MS_consensus_creation](../MS_consensus_creation/README.md).


## Important Notes
- To change parameter within the analysis edit the NASE.ini or decoy_database.ini . The details of each parameter are explained on the linked [OpenMS](../OpenMS/)
