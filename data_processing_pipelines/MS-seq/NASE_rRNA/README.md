# NASE rRNA
This folder contains the rRNA-specific wrapper for the Human RNome MS-seq search workflow. It reproduces the relevant OpenMS DecoyDatabase to NucleicAcidSearchEngine
steps on an `.mzML` file and an RNA FASTA reference, then writes a fragments containing `bedRMod` output file suitable for downstream consensus creation.

The script in this folder is the rRNA analogue of the tRNA workflow and is designed to work with the bundled configuration files in this directory.


## Installation
A conda package will be available soon. For now install openms via following .deb file https://cloud.samwein.com/s/EMizZJcsiZKsKfE OpenMS-3.6.0-pre-add-bedrmod-2026-04-14-Debian-Linux-x86_64.deb. The package was only tested on Linux-Debian version X and might not work on other systems.

**Install the RNOME specific OpenMS version:**
    ```bash
    sudo apt install OpenMS-3.6.0-pre-add-bedrmod-2026-04-14-Debian-Linux-x86_64.deb
    ```

**Create a python enviroment**
1. Ensure you have Conda (or Miniconda/Mamba) installed.
2. Create the environment by running:
   ```bash
   conda create -n rnome_openms python=3.10
   ```
3. Activate the environment:
   ```bash
   conda activate rnome_openms
   ``` 
4. *optional donwload package for raw file conversion
    ```bash
    conda install -c bioconda -c conda-forge proteowizard
    ``

## Input data
The workflow expects:

- an input `.mzML` file containing your MS-seq data
- a FASTA sequence file for the target RNA species 
- configuration files in this folder:
  - `decoy_database.ini`
  - `NASE.ini`
  - `ChEBI_ID_RNA_mods_compatible.csv`

The default rRNA refrecne in this folder is:

- `hs_rRNAs.fa`

### Usage

**In case you have a non .mzML raw data file**
```bash
msconvert file_name.X --mzML --outfile file_name.mzML
```

**main pipeline**
```bash
python run_tRNA_analysis.py \
    <input.mzML> \
    <input.fasta> \
    --precursor-tolerance 10 \
    --product-tolerance 20 \
    --output-dir ./results \
    --decoy-ini ./decoy_database.ini \
    --nase-ini ./NASE.ini \
    --chebi-mapping ./ChEBI_ID_RNA_mods_compatible.csv
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


## Importent Notes

- The percursor and product tolerance values used in the RNome Paper are listed in the tolerances.ods file as an example
- To change paramter within the analysis edit the NASE.ini or decoy_database.ini . The details of each parameter are explained on the linked [OpenMS](../OpenMS/)
