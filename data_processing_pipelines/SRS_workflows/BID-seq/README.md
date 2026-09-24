# Processing of BID-seq Short-Read Sequencing Data (He_pU_input / He_pU_treat_rep1 / He_pU_treat_rep2)

End-to-end workflow for the pseudouridine (&Psi;) RNome samples of the Human RNome Project,
from raw FASTQ files to a bedRMod file.

---

# Data Download Instructions

1. Open the project's data repository and follow the instructions under "Raw data".
2. Download the paired-end Illumina NovaSeq X FASTQ files listed below into a dedicated
   directory on your machine or HPC cluster (referred to as `${rawdir}` below).

| HRP sample name    | Library      | Raw FASTQ files                                                                     |
| ------------------ | ------------ | ----------------------------------------------------------------------------------- |
| `He_pU_input`      | untreated    | `CHe-MA-42s-YSL-04_S4_L002_R1_001.fastq.gz`, `CHe-MA-42s-YSL-04_S4_L002_R2_001.fastq.gz` |
| `He_pU_treat_rep1` | BID-treated  | `CHe-MA-42s-YSL-07_S7_L002_R1_001.fastq.gz`, `CHe-MA-42s-YSL-07_S7_L002_R2_001.fastq.gz` |
| `He_pU_treat_rep2` | BID-treated  | `CHe-MA-42s-YSL-08_S8_L002_R1_001.fastq.gz`, `CHe-MA-42s-YSL-08_S8_L002_R2_001.fastq.gz` |

Input RNA: poly(A)-selected, fragmented human mRNA. Platform: Illumina NovaSeq X, paired-end.

---

# Requirements & Setup

The complete analysis runs inside a published container, so no conda environment is required
for the main pipeline.

| Software                    | Version                      | Purpose                                     |
| --------------------------- | ---------------------------- | ------------------------------------------- |
| apptainer (or singularity)  | >= 1.1                       | container runtime                           |
| `docker://y9ch/bidseq`      | pipeline **v2.0**            | the entire BID-seq workflow                 |
| STAR                        | 2.7.7a                       | building the genome index (reference prep)  |
| samtools                    | >= 1.17                      | extracting the rRNA reference               |
| python                      | >= 3.8                       | `convert_to_bedRmod.py`                     |

The container image is built from <https://github.com/y9c/pseudoU-BIDseq>
(source revision used here: commit `750f995`, 2024-05-08, `VERSION` = `v2.0`) and bundles
fastp, cutadapt, bowtie2, STAR, samtools, bedtools, UMICollapse and snakemake 7.18.

```sh
# on a node with internet access
module load apptainer          # if not available by default
apptainer pull bidseq.sif docker://y9ch/bidseq

# record the exact image that was used
apptainer inspect bidseq.sif
```

---

# Input

- Paired-end FASTQ files (`R1` + `R2`) as downloaded; no pre-trimming is required.
  Adapter removal, UMI handling, read joining and deduplication are all performed by the
  pipeline.
- The libraries carry a 5-nt random UMI at the 5' end of read 2
  (pipeline setting `barcode: '-NNNNN'`, which is the default).
- The libraries are reverse-stranded: read 1 is antisense to the transcript
  (`forward_stranded: false` in the configuration file).

---

# Workflow

## 1. Prepare the references

Two references are needed: the **genome** (for transcriptome-wide &Psi; calling) and a small
**`genes`** reference (for the RNome-specific rRNA analysis, which is mapped before the genome
and therefore recovers the multi-copy 45S rDNA reads that are lost in genome mapping).

```sh
refdir='your_reference_directory'
```

### 1a. Genome reference and STAR index

Ensembl GRCh38 primary assembly with the release-110 annotation:

```sh
curl -O https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
curl -O https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > ${refdir}/GRCh38.fa
gunzip -c Homo_sapiens.GRCh38.110.gtf.gz               > ${refdir}/GRCh38.release110.gtf
samtools faidx ${refdir}/GRCh38.fa

STAR --runMode genomeGenerate --runThreadN 24 \
     --genomeDir ${refdir}/star/GRCh38.release110 \
     --genomeFastaFiles ${refdir}/GRCh38.fa \
     --sjdbGTFfile ${refdir}/GRCh38.release110.gtf \
     --genomeSAindexNbases 14 \
     --limitGenomeGenerateRAM 55000000000
```

(STAR 2.7.7a was used; `sjdbOverhang` is left at its default of 100.)

### 1b. RNome `genes` reference (human 45S pre-rRNA)

The `genes` reference used for these samples is a single sequence: the human 45S pre-rRNA
unit **RNA45SN1** (GeneID 106631777), located on GRCh38 chromosome 21 at
**8,433,222-8,446,572 (+ strand)**.

```sh
samtools faidx ${refdir}/GRCh38.fa 21:8433222-8446572 \
  | sed '1s/.*/>NC_000021.9:8433222-8446572 RNA45SN1 [organism=Homo sapiens] [GeneID=106631777]/' \
  > ${refdir}/rRNA.fa
samtools faidx ${refdir}/rRNA.fa
```

The FASTA header matters: its first field becomes the `chr` value in
`filter_sites/genes.tsv.gz`, and `convert_to_bedRmod.py` uses it to lift the rRNA-local
coordinates back onto chromosome 21 (see step 4).

The pipeline builds its own bowtie2 index for this reference automatically; no manual index
build is needed. No `contamination` reference was used for these samples.

## 2. Write the configuration file

Save the following as `data.yaml` (a copy is provided next to this README). Adjust the paths.

```yaml
reference:
  genes:
    fa: /path/to/ref/rRNA.fa
  genome:
    fa: /path/to/ref/GRCh38.fa
    star: /path/to/ref/star/GRCh38.release110

samples:
  He_pU_input:
    data:
      - R1: /path/to/raw/CHe-MA-42s-YSL-04_S4_L002_R1_001.fastq.gz
        R2: /path/to/raw/CHe-MA-42s-YSL-04_S4_L002_R2_001.fastq.gz
    group: HRP_pU
    treated: false
  He_pU_treat_rep1:
    data:
      - R1: /path/to/raw/CHe-MA-42s-YSL-07_S7_L002_R1_001.fastq.gz
        R2: /path/to/raw/CHe-MA-42s-YSL-07_S7_L002_R2_001.fastq.gz
    group: HRP_pU
    treated: true
  He_pU_treat_rep2:
    data:
      - R1: /path/to/raw/CHe-MA-42s-YSL-08_S8_L002_R1_001.fastq.gz
        R2: /path/to/raw/CHe-MA-42s-YSL-08_S8_L002_R2_001.fastq.gz
    group: HRP_pU
    treated: true

forward_stranded: false
```

All three libraries share one `group`, because the two treated replicates and the untreated
control describe the same biological sample; the `group` is the unit over which sites are
pre-filtered and over which the &Psi; fraction is computed.

> Note on column names: in the original run the three samples were named `10ngut`,
> `10ngrep1` and `10ngrep2` with `group: 10ng`, which is why the shipped intermediate table
> `BIDSEQ_rRNA_rep1_2.tsv` carries the column prefixes `10ngut_*`, `10ngrep1_*`, `10ngrep2_*`
> and `10ng_ratio` / `10ng_fraction` / `10ng_passed`. The sample names above are the HRP
> names; only the column prefixes in the output change.

## 3. Run the whole analysis with one command

```sh
apptainer run -B /the/real/path bidseq.sif -c data.yaml -j 48
```

`-B` must bind the real (resolved) path of every directory holding the references, the FASTQ
files and the working directory; use `realpath ./` to find it on shared HPC file systems.
`-j` sets the number of parallel jobs/cores (48 was used).

A SLURM submission script for this exact run is provided as `run_bidseq.sbatch`.

### Processing and site-calling parameters

All parameters below are the pipeline defaults (`config.yaml` inside the container) and were
used unchanged for these samples:

```yaml
adapter:
  p7: AGATCGGAAGAGCACACGTC
  p5: AGTTCTACAGTCCGACGATC
barcode: '-NNNNN'          # 5-nt UMI at the 5' end of read 2
trim_p5: false
trim_polyA: false
greedy_mapping: false
speedy_mapping: false

cutoff:                    # pre-filtering of candidate sites
  min_match_prop: 0.8      # STAR outFilterMatchNminOverLread
  min_group_gap: 5         # >= 5 deletions summed over the group
  min_group_depth: 10      # >= 10x coverage summed over the group
  min_group_ratio: 0.01    # >= 1% deletion ratio in the group
  min_group_num: 1

group_filter:              # final site filter
  combine_group_input: true
  min_passed_group: 1
  min_treated_depth: 20
  min_input_depth: 20
  min_treated_gap: 5
  min_treated_ratio: 0.02
  min_treated_fraction: 0.02
  min_fold_ratio: 2        # treated deletion ratio >= 2x the input
  max_p_value: 0.0001
```

&Psi; stoichiometry is obtained by converting the observed deletion ratio with the calibration
curves shipped in the container (`/pipeline/calibration_curves.tsv`).

## 4. Output of the pipeline

Three folders are created in `workspace/`:

```
workspace/
├── align_bam/      # {sample}_genes.bam, {sample}_genome.bam (+ .bai)
├── report_reads/   # trimming / mapping / deduplication reports (MultiQC html)
└── filter_sites/
    ├── genes.tsv.gz    # sites on the 45S rRNA reference
    └── genome.tsv.gz   # sites on GRCh38
```

Both site tables share one format:

| Column               | Meaning                                                                    |
| -------------------- | -------------------------------------------------------------------------- |
| `chr`                | sequence name (`1`, `2`, ... for the genome; the FASTA header for `genes`)  |
| `pos`                | 1-based position of the U                                                  |
| `strand`             | `+` / `-`                                                                  |
| `{sample}_depth`     | read coverage of that library at the site                                  |
| `{sample}_gap`       | reads carrying the BID-seq deletion signature at the site                  |
| `{group}_ratio`      | deletion ratio of the treated libraries of the group                       |
| `{group}_fraction`   | calibrated &Psi; fraction (stoichiometry), 0-1                             |
| `{group}_passed`     | `1` if the site passes all `group_filter` thresholds                       |

## 5. Convert the results into bedRMod

`convert_to_bedRmod.py` (provided next to this README) turns either site table into a
bedRMod v2 file. It needs only a standard python >= 3.8 (no third-party packages).

```sh
# transcriptome-wide sites (GRCh38 coordinates)
python convert_to_bedRmod.py \
    --input  workspace/filter_sites/genome.tsv.gz \
    --output He_pU_mRNA.bedrmod \
    --treated-samples He_pU_treat_rep1,He_pU_treat_rep2 \
    --group HRP_pU --passed-only

# rRNA sites, lifted from the 45S reference onto chromosome 21
python convert_to_bedRmod.py \
    --input  workspace/filter_sites/genes.tsv.gz \
    --output He_pU_rRNA.bedrmod \
    --treated-samples He_pU_treat_rep1,He_pU_treat_rep2 \
    --group HRP_pU --passed-only \
    --lift 'NC_000021.9:8433222-8446572=21:8433222:+'
```

`--lift NAME=CHROM:START:STRAND` maps a local 1-based position `p` on the reference sequence
`NAME` to the genomic 1-based position `START + p - 1`. Use
`--chrom-style ucsc` to emit `chr21` instead of `21`.

Conventions written by the script:

- `chromStart` is 0-based (`pos - 1`), `chromEnd = chromStart + 1`;
- `name` is the modification short name `Y` (`--name-field id` writes `17802` instead);
- `score` is `min(1000, coverage)`;
- `coverage` is the summed depth of the treated libraries;
- `frequency` is the calibrated &Psi; fraction as a **percentage** (0-100, two decimals).

Both files can be concatenated after conversion (keep one header) if a single bedRMod file
per sample is required.

---

# Output

- `He_pU_mRNA.bedrmod` - &Psi; sites transcriptome-wide, GRCh38 coordinates.
- `He_pU_rRNA.bedrmod` - &Psi; sites on the 45S pre-rRNA, reported on chromosome 21.

Header of the produced files:

```
#fileformat=bedRModv2
#organism=9606
#modification_type=RNA
#modification_names=17802:Y:U
#assembly=GRCh38
#annotation_source=Ensembl
#annotation_version=110
#sequencing_platform=Illumina NovaSeq X
#bioinformatics_workflow=pseudoU-BIDseq v2.0 (docker://y9ch/bidseq)
```
