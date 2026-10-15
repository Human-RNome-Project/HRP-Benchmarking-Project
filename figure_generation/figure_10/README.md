# Figure 10 b,c: short-read sequencing (SRS) modification maps

This folder contains the plotting scripts used to generate Figure 10b and 10c in the manuscript.

## Overview

The folder contains two main panel generators:

- `figure10b.py` — assembles the rRNA lollipop-style maps for `5.8S`, `18S`, and `28S` rRNA.
- `figure10c.py` — builds the tRNA heatmap using Sprinzl coordinate alignment.

The rRNA panel rendering logic is split into three dedicated plotting helpers:

- `plot_18.py` — 18S lollipop map
- `plot28.py` — 28S lollipop map
- `plot_5_8.py` — 5.8S lollipop map

## Contents

- `figure10b.py`
  - Main entry script for the rRNA panels.
  - Loads the GM12878 SRS consensus BED file and compares it against a reference file.
  - Calls the `plot_18`, `plot_28`, and `plot_5_8` functions to render the three rRNA maps.

- `figure10c.py`
  - Generates the tRNA Illumina heatmap.
  - Reads the GM12878 SRS consensus BED file, normalizes tRNA naming, aligns positions using Sprinzl coordinates, and builds a per-position modification matrix.
  - Writes an output SVG named `tRNA_Illumina_heatmap.svg`.

- `plot_18.py`
  - Plotting helper for the `18S` lollipop map.
  - Colors and orders modification calls by modification type and overlap status.

- `plot28.py`
  - Plotting helper for the `28S` lollipop map.
  - Uses the same visual logic as the 18S helper, with a different subplot size/layout.

- `plot_5_8.py`
  - Plotting helper for the `5.8S` lollipop map.
  - Uses the same modification-level lollipop style as the other rRNA panels.

- `tRNA_sprinzl.xlsx`
  - Reference workbook used to map tRNA positions into Sprinzl coordinates for the tRNA heatmap.

- `H.sapiens_rRNA_reference.bed`
    - rRNA modification reference from literature, Taoka et al. 2018, adapted by Yuri Motorin.

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

### rRNA 10b
- a reference rRNA BED file for comparison (we used Taoka et al. 2018, adapted by Yuri Motorin)
- a rRNA consensus BED file containing the consensus modification calls

### tRNA 10c
- a tRNA consensus BED file containing the consensus modification calls
- the `tRNA_sprinzl.xlsx` coordinate table for tRNA alignment

## Typical Workflow

1. Prepare the SRS consensus BED data for rRNA and tRNA.
2. Adapt the main function in each script to point to your data files.
3. Run the panel entry script:
   - `python figure10b.py` for the rRNA lollipop maps
   - `python figure10c.py` for the tRNA heatmap
4. Review the generated SVG outputs in the working directory.


