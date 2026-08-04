# Figure 6: MS-seq modification maps

This folder contains the plotting scripts used to generate Figure 6 in the manuscript.

## Overview

Figure 6 combines three panels:

- `figure_6a.py` — plots the 18S rRNA MS-seq modification map.
- `figure_6b.py` — plots the 28S rRNA MS-seq modification map.
- `figure_6c.py` — plots the tRNA modification heatmap for the 47 isoacceptor families.

A shared helper module, `figure_6_helper_functions.py`, provides the common data-loading, modification-name normalization, fragment extraction, and color-mapping logic used by the panel scripts.

## Contents

- `figure_6a.py`
  - Generates the 18S rRNA panel.
  - Uses a detailed consensus rRNA BED file plus a Taoka reference BED file to identify and visualize overlapping vs non-overlapping modifications.
  - Writes an SVG output named `rRNA_18S_modifications.svg`.

- `figure_6b.py`
  - Generates the 28S rRNA panel.
  - Mirrors the same plotting logic as the 18S script, but for the 28S region.
  - Writes an SVG output named `rRNA_28S_modifications.svg`.

- `figure_6c.py`
  - Builds the tRNA heatmap showing per-position modification coverage and modification identity across tRNA structural regions.
  - Uses the Sprinzl coding coordinates to align positional annotations across all tRNAs.
  - Produces a large tRNA heatmap figure, suitable for publication-style visualization.

- `figure_6_helper_functions.py`
  - Shared helper code for:
    - chromosomal name normalization
    - modification code mapping
    - BEDRmod loading and filtering
    - fragment extraction and reference overlap logic
    - consistent modification color definitions

- `tRNA_sprinzl.xlsx`
  - Reference workbook used for tRNA Sprinzl coordinate alignment in the tRNA heatmap plot.

- `H.sapiens_rRNA_reference.bed`
    - rRNA modification reference from literature, Taoka et al. 2018, adapted by Yuri Motorin

- `environment.yml`
    - a .yml file containing the used conda environment.



## Enviroment setup
The scripts require `pandas`, `numpy`,`openpyxl` and `matplotlib`.

### Using Conda
1. Ensure you have Conda (or Miniconda/Mamba) installed.
2. Create the environment by running:
   ```bash
   conda env create -f environment.yml
   ```
3. Activate the environment:
   ```bash
   conda activate MS-seq_consensus_creation
   ```

## Expected Inputs

The file paths in `main()` blocks need to be adapted to your local system before usage. In practice, they expect:

### rRNA 6a,b
- a reference rRNA BED file for comparison (we used Taoka et al. 2018, adapted by Yuri Motorin)
- a rRNA consensus BED file containing the consensus modification calls (created via the [MS consensus creation script](../../data_processing_pipelines/MS_consensus_creation/)) 
- raw rRNA `.bed` files from [OpenMS](../../data_processing_pipelines/OpenMS/)

### tRNA 6c
- a tRNA consensus BED file containing the consensus modification calls (created via the [MS consensus creation script](../../data_processing_pipelines/MS_consensus_creation/))
- raw tRNA `.bed` files from [OpenMS](../../data_processing_pipelines/OpenMS/)
- a tRNA BEDRmod dataset and matching Sprinzl coordinate table


## Typical Workflow

1. Prepare the input BEDRmod data files.
2. Adapt the main function in each script to point to your data files.
3. Run the appropriate panel script:
   - `python figure_6a.py`
   - `python figure_6b.py`
   - `python figure_6c.py`
4. Inspect the generated SVG/PDF-style output in the current working directory.


