# Human RNome Project — Data Processing Pipelines

This directory contains the computational data processing pipelines used by the **Human RNome Project (HRP)** to process raw sequencing and mass spectrometry datasets, perform quality control and statistical validation, build consensus modification references, and prepare genome browser visualization hubs.

---

## Overall Pipeline Architecture

The overall data flow, from the raw data inputs (`DOE Data Explorer doi:10.25585/DOE-HRP/3377574`) through the processing pipelines in this folder to final figure generation, is illustrated below (corresponding to `HRP-benchmarking-project/flowchart.png`):

```
graph TD
    Raw[Raw Data Inputs<br/>DOE Data Explorer] --> LRS_Proc[1. LRS_modification_pipeline]
    Raw --> MS-seq[2. MS-seq]
    Raw --> SRS_workflows[SRS_workflows]

    LRS_Proc -->|Primary Basecalling & Validation| LRS_Bed[LRS Platform BED File]
    MS-seq -->|Peak Identification & Spectra| MS_Cons[3. MS_consensus_creation]
    MS_Cons -->|Filtering & Sample Merging| MS_Bed[MS Platform BED File]

    LRS_Bed --> RDD[4. RNA_modification_estimation]
    RDD -->|Error Rates & RDD Metrics| Fig13[Figure 13 Plotting<br/>figure_generation/figure_13_RNA_modification_estimation]

    LRS_Bed --> GlobalCons[5. consensus-draft-sequence]
    SRS_Bed --> GlobalCons
    MS_Bed --> GlobalCons

    GlobalCons -->|Multi-platform Jaccard Grid Search| ConsFile[Consensus Draft Sequence File<br/>consensus_draft_sequence.bed]
    ConsFile --> Fig12[Figure 12 Consensus Plotting<br/>figure_generation/figure_12_consensus]

    LRS_Bed -.-> UCSC[6. HRP_UCSC_hub]
    SRS_Bed -.-> UCSC
    MS_Bed -.-> UCSC
    ConsFile --> UCSC
    UCSC --> UCSCHub[UCSC Genome Browser Track Hub]
```

---

## Directory Overview

The processing workflows in this directory are structured into three platform-specific pipeline components that need to be implemented in no particular order:

### 1\. [`LRS_modification_pipeline/`](http://LRS_modification_pipeline/) (Primary Long-Read Sequencing Data Processing)

- **Purpose:** Processes raw Oxford Nanopore Direct RNA Sequencing (DRS) POD5 files into statistically validated RNA modification calls across polyA RNA, rRNA, and tRNA biotypes.
- **Key Modules:**
  - **`wdl-pipelines-main/`:** Upstream WDL workflows for Dorado basecalling, Minimap2 alignment, NanoComp QC, and Modkit modification pileup.
  - **`format_conversion_and_merging/`:** Standardizes raw Modkit outputs from **`wdl-pipelines-main/`** into `bedRMod` format, applies coordinate adapter offsets, and handles transcriptome liftovers to create a DRS platform RNA reference.
  - **`wf-LRS-bedrmod-native-vs-ivt-main/`:** Contains the Snakemake workflow performing Fisher's exact tests comparing native vs. in-vitro transcribed (IVT) control samples with Benjamini-Hochberg FDR correction.
- **Flowchart Integration:** Input: `Raw POD5 Data` → Output: `LRS Platform BED File` (feeds into `consensus-draft-sequence`, `RNA_modification_estimation`, `HRP_UCSC_hub`, and Figures 8 & 9).

---

### 2\. [`MS-seq/`](https://github.com/Human-RNome-Project/HRP-Benchmarking-Project/tree/main/data_processing_pipelines/MS-seq) (Primary Mass Spectrometry Data Processing)

- **Purpose:** Contains toolset for processing liquid chromatography-mass spectrometry (LC-MS/MS) nucleoside and oligonucleotide datasets. Also merges, filters, and standardizes sample-level Mass Spectrometry (OpenMS) bedRmod files into a single unified MS platform reference.
- **Key Modules:**
  - **`NASE_rRNA/`:** contains the rRNA-specific wrapper for the Human RNome MS-seq search workflow to write `bedRmod` output files suitable for downstream `MS_consensus_creation` step.
  - **`NASE_tRNA/`:** contains the tRNA-specific wrapper for the Human RNome MS-seq search workflow to write `bedRmod` output files suitable for downstream `MS_consensus_creation` step.
  - **`OpenMS/` :** Open-source C++ framework and toolset for processing liquid chromatography-mass spectrometry (LC-MS/MS) nucleoside and oligonucleotide datasets.
  - **`MS_consensus_creation/`:** Upstream WDL workflows for Dorado basecalling, Minimap2 alignment, NanoComp QC, and Modkit modification pileup.
  - **`docker/`:** Provides a docker entry point bash script to utilize the Docker image and MS-seq functions.
- **Flowchart Integration:** Input: Sample OpenMS BEDs → Output: `MS Platform BED File` (feeds into `consensus-draft-sequence`, `HRP_UCSC_hub`, and Figures 3 & 6).

---

Once all platform-specific processing pipelines (SRS\_workflows, LRS\_modification\_pipeline, MS-seq) are complete, we will create a consensus RNA sequence.

### 1\. [`consensus-draft-sequence/`](http://consensus-draft-sequence/) (Global Consensus Draft Reference Integration)

- **Purpose:** Integrates orthogonal modification calls across Short-Read Sequencing (SRS), Long-Read Sequencing (LRS), and Mass Spectrometry (MS) platforms to build the Human RNome Project benchmark consensus draft reference.
- **Key Modules:**
  - **`scripts/harmonize_massspec.py`:** Resolves ambiguous MassSpec modification calls (`mxA/mxC/mxG/mxU`) against curated reference BEDs and high-confidence sequencing datasets.
  - **Grid Search & Optimization:** Quantile parameter range generation (`scripts/generate_parameter_ranges.py`) and multi-way grid search optimization (`scripts/run_rrna_grid_from_ranges.py`, `scripts/run_3way_*.py`, `scripts/run_trna_grid_from_ranges.py`) to maximize Jaccard similarity across platforms.
  - **Classification & Tiering:** Filtering (`run_filtering.py`) and creation of tiered modification lists (`scripts/create_tiered_mod_lists.py`) and the unified `consensus_draft_sequence.bed`.
- **Flowchart Integration:** Inputs: `LRS Platform BED`, `SRS Platform BED`, `MS Platform BED` → Output: `Consensus Draft Sequence File` (`consensus_draft_sequence.bed`) (feeds into **Figure 12** plotting in `figure_generation/figure_12_consensus/` and `HRP_UCSC_hub`).

---

The following two pipelines correspond to estimating RNA modifications using sequencing and building a UCSC Genome Browser track to display the consensus RNA sequence.

### 1\. [`RNA_modification_estimation/`](http://RNA_modification_estimation/) (Secondary Error Rate & RDD Analysis)

- **Purpose:** Analyzes sequencing error profiles, substitution spectra, and candidate RNA-DNA Differences (RDDs) from long-read alignments.
- **Key Modules:**
  - Upstream per-chromosome `pysamstats --type variation` count table generation (`samtools view -F 2304`).
  - **`scripts/generation/`:** Upstream calculation scripts for per-read total error rates (`calculate_read_error_rates_parallel.py`), coverage-binned mismatch rates (`calculate_coverage_error_rate.py`), modkit stratified mismatch rates (`calculate_native_mismatch_rates_memeff_fast.py`), and genome-wide RDD JSON summaries (`rdd_genomewide.py`).
- **Flowchart Integration:** Input: `LRS Platform BED` & aligned BAMs → Output: Precalculated error & RDD summary tables (feeds directly into **Figure 13** plotting in `figure_generation/figure_13_RNA_modification_estimation/`).

---

### 2\. [`HRP_UCSC_hub/`](http://HRP_UCSC_hub/) (Genome Browser Track Hub Generation)

- **Purpose:** Converts modification call datasets into bigBed tracks and track hub configurations for interactive visualization on the UCSC Genome Browser.
- **Key Modules:**
  - **`scripts/pipeline.py`:** Main orchestration script splitting BED files by RNA biotype, building bigBed files, and generating track hub configurations.
  - **`scripts/make_decorator.py`:** Generates modification decorator tracks overlaid on underlying RNA sequence tracks.
  - **`scripts/BedPyLift.py`:** Utility for transcriptome-to-genome coordinate liftovers.
- **Flowchart Integration:** Inputs: `LRS Platform BED`, `SRS Platform BED`, `MS Platform BED`, `Consensus Draft Sequence File` → Output: Ready-to-load UCSC Track Hub (`ucsc_hub/`) for web browsing and Figure 14\.

---

## Mapping Pipelines to Figure Generation

| Data Processing Directory | Key Output File(s) | Primary Downstream Figure Target |
| :---- | :---- | :---- |
| **`LRS_modification_pipeline/`** | `LRS Platform BED File` | `figure_generation/figure_08`, `figure_09` |
| **`MS-seq/`** | `MS Platform BED File` | `figure_generation/figure_03`, `figure_06` |
| **`RNA_modification_estimation/`** | Count tables & RDD summaries | `figure_generation/figure_13_RNA_modification_estimation/` |
| **`consensus-draft-sequence/`** | `consensus_draft_sequence.bed`, `tiered_*.tsv` | `figure_generation/figure_12_consensus/` |
| **`HRP_UCSC_hub/`** | UCSC Track Hub (`ucsc_hub/`) | Genome Browser / Figure 14 |
