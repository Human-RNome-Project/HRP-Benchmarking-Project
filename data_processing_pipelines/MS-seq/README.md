# MS-seq Analysis

This directory contains the pipeline modules used to process MS-seq data and identify post-transcriptional RNA modifications across the RNA biotypes rRNA and tRNA.

The pipeline is organized into two distinct layers, the Open-MS spectra to a sequence mapping and the Consensus creation pipeline, combining aligned fragments and different samples to a consensus modification map.

---

## Overall Pipeline Architecture

The workflow progresses through the folders in this directory as follows:

```mermaid
graph TD
    A[Raw Data (.mzML)] --> B[1. OpenMS NucleicacidSearchEngine (NASE)/ NASE_rRNA or NASE_tRNA]
    B -->|Fragment spectra generation + Spectra Matching + Scoring| C[Aligned Fragments .bedrmod]
    C --> D[2. MS_consensus_creation]
    D -->|Fragment Overlap + Sample Overlap + Conditional Filtering| E[Consensus Modifiaction File .bed + Stastics Consensus Modifiaction File .bed ]
```

---


## Directory Overview

The folder is divided into four main components:

### 1. [`NASE_rRNA/`](NASE_rRNA/) (rRNA-pipeline)
- **Purpose:** rRNA specific generation of digested fragments and their spectra, matching of these spectras with the sample, scoring of the results 
- **Key Features:**
    - Pipeline excecution script (see details in [README.md](NASE_rRNA/README.md)).
    - Pipeline Specifcations for generation of decoy database and general workflow

### 1. [`NASE_tRNA/`](NASE_tRNA/) (tRNA-pipeline)
- **Purpose:** tRNA specific generation of digested fragments and their spectra, matching of these spectras with the sample, scoring of the results 
- **Key Features:**
    - Pipeline excecution script (see details in [README.md](NASE_rRNA/README.md)).
    - Pipeline Specifcations for generation of decoy database and general workflow


### 3. [`MS-seq Consensus creation/`](MS_consensus_creation/) (Sample Merging/Filtering pipeline)
- **Purpose:** Filtering and merging of different aligned fragment files to a consensus modification map.
- **Key Features:**
  - Consensus_creation_main python script. Initiate the pipeline (see details in [README.md](MS_consensus_creation/README.md)).

### 4. [`OpenMS/`](OpenMS/) (Main MS software bundle)
- **Purpose:** Links to the OpenMS github page. Further Details for devlopers and scientists to adapt our pipeline and use additional tools.

---

## Getting Started

1. **rRNA/tRNA Processing:** Follow the instructions in the [`NASE_rRNA/README.md`](NASE_rRNA/README.md) or [`NASE_tRNA/README.md`](NASE_tRNA/README.md) to get the aligned modified fragments.
2. **Filtering and Consensus Modifications** Follow the instructions in the [`MS_consensus_creation/README.md`](MS_consensus_creation/README.md) to extract your modifications.




