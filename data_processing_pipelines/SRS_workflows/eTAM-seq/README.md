# eTAM-seq — Data Processing Pipeline

This directory contains the computational data processing pipelines used by the **Human RNome Project (HRP)** to process evolved TadA-assisted N6-methyladenosine sequencing (eTAM-seq).

---
## Directory Overview

The processing workflows in this directory are structured into two specialized pipeline components. In the following order, run:

### 1. [`PreProcess/`]
- **Purpose:** Processes raw eTAM-seq FASTQ files into annotated pileup tables. The output produced here will be used in the /eTAM-seq/DownstreamProcess/pipeline.R.

---
### 2. [`DownstreamProcess/`] 
- **Purpose:** This directory contains the R pipeline used to quantify m6A sites from the eTAM-seq annotated tables from Step 1. This directory also compares methylation levels between mRNA samples and their in vitro transcribed control.

---
