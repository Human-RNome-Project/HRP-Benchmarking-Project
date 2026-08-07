# Figure 3: MS Quantification

This folder contains the plotting scripts used to generate Figure 3 in the manuscript.

## Overview

Figure 3 is assembled from two related analysis workflows:

- `figure3a-f.py` generates the MS quantification panels for Figure 3a–f. These panels compare modification abundance across rRNA fractions (5S, 5.8S, 18S, 28S, tRNA and total RNA) using the intermediate Excel tables provided in this folder.
- `figure_3g.py` generates the cap-structure bar plot shown in Figure 3g.

## Contents

- `figure3a-f.py`
  - Loads the intermediate Excel exports for HRP-C-003, HRP-C-004/005/006/007/008/0010, HRP-C-021/022 and the combined mean/std workbook.
  - Normalizes modification names with `translation_dicts.py`.
  - Produces for each analyzed RNA type a SVG barplot output.

- `figure_3g.py`
  - Loads `HRP_C_009_cap_structures.xlsx`.
  - Converts the abundances to the requested scale and plots the cap-structure composition as a broken-axis bar plot.
  - Writes `cap_modifications_barplot.svg`.

- `general_plot.py`
  - Shared plotting utilities for the broken-axis visualizations used in figure 3a–f.

- `translation_dicts.py`
  - Translation table used to standardize modification names across the input sheets.

- `environment.yml`
  - Conda environment specification for the required plotting dependencies.

## Environment setup

The scripts require Python 3.11 and the following libraries:

- `pandas`
- `numpy`
- `openpyxl`
- `matplotlib`

### Using Conda

1. Ensure Conda or Mamba is installed.
2. Create the environment from the provided file:
   ```bash
   conda env create -f environment.yml
   ```
3. Activate the environment:
   ```bash
   conda activate figure_env
   ```

## Expected inputs

The scripts expect the Excel files in this directory to be available with the names referenced inside the scripts. In practice, these files are used:

- `HRP-C-003_intermediate_results.xlsx`
- `HRP-C-004_005_006_007_008_0010_intermediate_results.xlsx`
- `HRP-C-021_022_intermediate_results.xlsx`
- `HRP-C-mean_std_combined.xlsx`
- `HRP_C_009_cap_structures.xlsx`

If these files are moved to a different location, the paths at the top of each script should be updated accordingly.

## Typical workflow

1. Open the folder and confirm that the expected Excel files are present.
2. Activate the conda environment described above.
3. Run the figure generation scripts from this directory:
   ```bash
   python figure3a-f.py
   python figure_3g.py
   ```
4. Inspect the generated SVG files in the current working directory.


