# eTAM-seq Pre-processing Pipeline

This directory contains the pre-processing chain that turns raw paired-end eTAM-seq FASTQ files
into the per-site **annotated pileup tables** (`<sample>.annot.txt`) consumed downstream by the
R pipelines (`/pipeline.R`).

The chain does three things: it cleans and stitches the read pairs (adapters, UMI, overlap
merging), it maps the stitched reads with **hisat-3n** in A→G three-letter mode and deduplicates
them by UMI, and it finally counts, at every adenosine of the transcriptome, the reads that still
read **A** (protected, i.e. methylated) versus those converted to **G**, adding the genomic
annotation (gene, feature kind, codon, sequence context, splice-site distance) taken from the
RefSeq GenBank flat files.

## Directory Structure

Repository side (`/PreProcess/`):

- `run`                                # Master driver: loops over the samples, trims, stitches, calls the two other scripts
- `2_map_and_count_stitched_reads.sh`  # rRNA + genome mapping with hisat-3n, dedup, conversion filter, pileup2var
- `annotPileupRev`                     # Perl: strand-aware mpileup A/G counting + GenBank annotation
- `README.md`                          # Contains project documentation

Data side (the working directory where `run` is launched, one level **below** this folder):

- `fq/`                                # Input FASTQ files, `<sample>_L?_R[12]_???.fastq.gz`
- `outstitch/<sample>/`                # BAMs, stats, pileup2var tables (`$hisat3n_out_dir`)
- `<sample>.annot.txt`                 # Final annotated table, written in the working directory
- `<sample>.a3a5trim.umi.round[12].log`, `<sample>.a3a5trim.umi.log`  # cutadapt / umi_tools logs

Reference side (`$REF_DIR`, **parameter**; default `/home/h/OneDrive/Unil/data/b38`):

- `b38.fa.gz`                          # Genome sequence, bgzipped (mpileup / annotPileupRev)  ✔ present
- `GB/*.gb.gz`                         # RefSeq GenBank flat files, one per chromosome         ✔ present
- `b38.fa`                             # Genome sequence, plain (samtools stats, pileup2var)   ✘ to recreate
- `HS_rep/b38*`                        # hisat-3n genome index                                 ✘ to rebuild
- `HS_rep/hs_rrna*`, `hs_rrna.fa`      # hisat-3n rRNA index and its FASTA                     ✘ to rebuild

The reference folder : see *Reference data: where to get it* below for
what the local copy contains, how to obtain the missing files, and the variables that point to
them.

### Expected layout at run time

The three files of this folder are copies of the scripts run on the server, where they sit on
**two levels**: `run` and `2_map_and_count_stitched_reads.sh` in the dataset folder (together
with the input `fq/`), and `annotPileupRev` one level above, shared by the other datasets.

```
<parent>/
├── annotPileupRev
└── <dataset>/
    ├── run
    ├── 2_map_and_count_stitched_reads.sh
    └── fq/
```

Both layouts work: `run` uses `./annotPileupRev` when it finds an executable copy in the
dataset folder, and falls back to `../annotPileupRev` otherwise (`ANNOT_PL` overrides both).

## Scripts/Workflows

| File | Purpose |
|------|---------|
| `run` | Sample loop. Two rounds of `cutadapt`, UMI extraction (`umi_tools extract`), read stitching (`NGmerge`), then calls `2_map_and_count_stitched_reads.sh` and `annotPileupRev`. Skips a sample whose final filtered BAM already exists, and `break`s the whole loop at the first error. |
| `2_map_and_count_stitched_reads.sh` | Takes one sample name as `$1`. Maps the stitched reads on the rRNA index then on the genome with `hisat-3n --base-change A,G`, keeps uniquely mapped reads (MAPQ ≥ 60), deduplicates on the UMI, filters reads by their A→G conversion rate, and produces the `pileup2var` conversion tables. Deletes the intermediate BAMs at the end. |
| `annotPileupRev` | Reads the GenBank flat files (stdin or `.gb` arguments), builds a per-position annotation of every A of every transcript, then runs two strand-restricted `samtools mpileup` passes on the final BAM and prints one line per covered A site with its A/G counts and annotation. |

## Requirements & Setup

### Environment setup

No Conda environment or Docker image is used: the tools must be on the `PATH` of the machine
running the pipeline.

**Versions used for the reference run** — reproduce these before comparing new results to the
existing tables:

| Tool | Version | Used by | Note |
|------|---------|---------|------|
| `cutadapt` | `4.6.dev2+g3407ac0` | `run` | Development build, not a released 4.6; `-j 0` = all cores |
| `umi_tools` | `1.1.4` | `run` (extract), `2_map_…sh` (dedup) | |
| `NGmerge` | `0.3` | `run` | Read stitching |
| `hisat-3n` | `hisat2-align-s 2.2.1-3n-0.0.3` | `2_map_…sh` | Default path `$HOME/src/github/hisat-3n/hisat-3n`, overridden by `HISAT3N` |
| `samtools` | `1.16.1-78-g4b247ef` | `2_map_…sh`, `annotPileupRev` | Git build between 1.16.1 and 1.17; ≥ 1.12 is needed for the `-e` filter expressions on the `Yf`/`Zf` tags |
| `pileup2var` | `1.1.0` | `2_map_…sh` | Must be on the `PATH` |
| `perl` | system perl | `annotPileupRev` | Core modules only (`strict`, `warnings`) |

Version check:

```bash
cutadapt --version
umi_tools --version
NGmerge --version
~/src/github/hisat-3n/hisat-3n --version
samtools --version
pileup2var --version
```

The three files must be executable (`chmod +x run 2_map_and_count_stitched_reads.sh annotPileupRev`).

### Resources

`NCPUS` defaults to 24, for both NGmerge and the mapping script: the chain is written
for a 24-core server. The genome mapping and the `umi_tools dedup` step are the long ones
(tens of minutes to a few hours per sample); `annotPileupRev` reloads the **whole genome FASTA
for every chromosome**, which adds a significant I/O cost. Each `2_map_…sh` call starts by
`rm -rf`-ing its own `outstitch/<sample>/` folder, so a re-run restarts that sample from scratch.
Disk: the reference folder needs ~1 GB for `b38.fa.gz`, ~3 GB once uncompressed, plus the
hisat-3n indexes (tens of GB with the repeat index).

## Input data

### Data Download Instructions

FASTQ files are not in this repository. You can download the FASTQ files from the RNome Data Explorer "Short-read Sequencing" folder under the "Raw Data" folder. The files needed for this pipeline correspond to the sample "HRP_B_002". After downloading, these FASTQ files are expected in `$FQ_DIR` (`fq/` by default) of the
working directory, with the Illumina-style naming used by the glob patterns of `run`. The
reference files come from `$REF_DIR` — where to download or rebuild each of them is detailed in
*Reference data: where to get it*.

| Input File | Description | File Type/Format | Script Usage |
|------------|-------------|------------------|--------------|
| `$FQ_DIR/<sample>_L?_R1_???.fastq.gz` | Read 1 of every lane of the sample | gzipped FASTQ | `run` (concatenated on the fly with `gzip -dc`) |
| `$FQ_DIR/<sample>_L?_R2_???.fastq.gz` | Read 2 (carries the 6 nt UMI at its 5′ end) | gzipped FASTQ | `run` |
| `$REF_DIR/HS_rep/hs_rrna*` | hisat-3n rRNA index | hisat-3n index | `2_map_…sh`, rRNA depletion pass |
| `$REF_DIR/hs_rrna.fa` | rRNA sequences | FASTA | `samtools stats`, `pileup2var` |
| `$REF_DIR/HS_rep/b38*` | hisat-3n genome index (GRCh38) | hisat-3n index | `2_map_…sh`, genome pass |
| `$REF_DIR/b38.fa` | Genome sequence, uncompressed | FASTA | `samtools stats`, `pileup2var` |
| `$REF_DIR/b38.fa.gz` | Genome sequence, bgzipped | FASTA.gz | `samtools mpileup` and sequence loading in `annotPileupRev` |
| `$REF_DIR/GB/*.gb.gz` | RefSeq GenBank records, one per chromosome | GenBank, gzipped | `annotPileupRev`, annotation source |

Both scripts check these files before doing anything and stop with an explicit message naming
the missing file and the current `REF_DIR`.

### Parameters to adapt

#### Paths: the `REF_DIR` parameter

Every reference path comes from a single variable, `REF_DIR`, which can be given three ways
(first one wins):

```bash
./run /path/to/b38                     # 1. first argument of run
REF_DIR=/path/to/b38 ./run             # 2. environment variable
./run                                  # 3. default: /home/h/OneDrive/Unil/data/b38
```

`run` exports `REF_DIR` (and `GENOME_INDEX_NAME`, `NCPUS`), so
`2_map_and_count_stitched_reads.sh` inherits it; called alone it also takes it as its second
argument: `./2_map_and_count_stitched_reads.sh <sample> /path/to/b38`.

Variables of `run` (all `${VAR:-default}`, so any of them can be exported before the call):

| Variable | Default | Meaning |
|----------|---------|---------|
| `REF_DIR` | `/home/h/OneDrive/Unil/data/b38` | Root of the reference data |
| `GENOME_FA_GZ` | `$REF_DIR/b38.fa.gz` | Bgzipped genome given to `annotPileupRev` |
| `GB_DIR` | `$REF_DIR/GB` | Folder of the `*.gb.gz` GenBank records |
| `GENOME_INDEX_NAME` | `b38` | Index name, also part of every output file name |
| `FQ_DIR` | `fq` | Input FASTQ folder |
| `SAMPLES` | the 8 `Rnome-*` samples | Space-separated sample list |
| `NCPUS` | `24` | Threads (NGmerge here, mapping/sorting downstream) |
| `MAP_SH` | `./2_map_and_count_stitched_reads.sh` | Mapping script |
| `ANNOT_PL` | `./annotPileupRev` if executable, else `../annotPileupRev` | Annotation script |

Variables of `2_map_and_count_stitched_reads.sh`:

| Variable | Default | Meaning |
|----------|---------|---------|
| `REF_DIR` (or 2nd argument) | `/home/h/OneDrive/Unil/data/b38` | Root of the reference data |
| `HISAT3N` | `$HOME/src/github/hisat-3n/hisat-3n` | hisat-3n executable |
| `REP_INDEX` / `REP_INDEX_NAME` | `$REF_DIR/HS_rep/hs_rrna` / `hs_rrna` | rRNA index |
| `GENOME_INDEX` / `GENOME_INDEX_NAME` | `$REF_DIR/HS_rep/b38` / `b38` | Genome index |
| `REP_FA` / `GENOME_FA` | `$REF_DIR/hs_rrna.fa` / `$REF_DIR/b38.fa` | Uncompressed FASTA |
| `FASTQ_DIR` | `.` | Where the `<sample>.stitched.fq.gz` are |
| `OUT_DIR` | `./outstitch` | Output folder |
| `NCPUS`, `A2G_PERCENT`, `COV`, `STRANDNESS` | `24`, `0.5`, `1`, `R` | See the table below |

`annotPileupRev` takes the genome FASTA as its second argument (what `run` does) or from the
`GENOME_FA_GZ` environment variable; and it stops if the FASTA is missing instead of failing later in the pileup.

#### Other parameters

In `run`:

| Parameter | Example | Meaning |
|-----------|---------|---------|
| `SAMPLES` | `Rnome-FTO1 … Rnome-mRNA-3` | List of the samples to process |
| `-g` / `-A` (round 1) | `GACGCTCTTCCGATCT` / `AGATCGGAAGAGCGTC` | 5′ adapter of R1 and 3′ adapter of R2 |
| `-a` / `-G` (round 2) | `NNNNNNAGATCGGAAGAGCACA` / `TGTGCTCTTCCGATCT` | 3′ adapter of R1 (preceded by the 6 nt UMI) and 5′ adapter of R2 |
| `-q 6`, `-m 0:46` | — | Quality trimming, and minimum length 0 for R1 / 46 for R2 |
| `--bc-pattern2="^(?P<umi_1>.{6}).*"` | 6 nt | UMI taken at the 5′ end of R2 |
| `NGmerge -d -e 20 -s -u 45 -n $NCPUS` | — | Dovetailed alignments allowed (≥ 20 nt overlap), shortest stitched read kept, input qualities up to 45 |
| GenBank / genome paths | `$GB_DIR/*.gb.gz`, `$GENOME_FA_GZ` | Annotation and sequence given to `annotPileupRev` |

In `2_map_and_count_stitched_reads.sh` (all the settings are in the `## main setting` block):

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `fastq_dir` | `.` | Where the `<sample>.stitched.fq.gz` files are |
| `hisat3n_rep_index` / `hisat3n_rep_index_name` | `$REF_DIR/HS_rep/hs_rrna` | rRNA index used for the depletion pass |
| `hisat3n_index` / `hisat3n_index_name` | `$REF_DIR/HS_rep/b38` | Genome index; the name is reused in every output file name |
| `hisat3n_out_dir` | `./outstitch` | Output folder, one sub-folder per sample |
| `rep_fa` / `genome_fa` | `$REF_DIR/hs_rrna.fa`, `$REF_DIR/b38.fa` | FASTA used by `samtools stats` and `pileup2var` |
| `ncpus` | 24 | Threads |
| `A2G_percent` | 0.5 | Minimum fraction of converted As (`Yf/(Yf+Zf)`) for a read to be kept |
| `cov` | 1 | Minimum coverage kept by `pileup2var` |
| `strandness` | `R` | Reverse-stranded library (stitched reads read the antisense strand) |

In `annotPileupRev`: the genome FASTA comes from the second command-line argument (as `run`
passes it) or from `GENOME_FA_GZ`. The minimum coverage to print a site (`nbA + nbG >= 10`) and
the 7 nt window of the `context` column are hardcoded.

### Reference data: where to get it

The folder used by default, `/data/b38`, is the local copy of the server's
`/data6/b38`. It currently holds **only the two files that cannot be regenerated quickly**:

| File | State | What it is |
|------|-------|------------|
| `b38.fa.gz` | present (≈ 970 MB, bgzip) | GRCh38 genome, UCSC-style names (`>chr1`, …) |
| `GB/NC_0000{01..24}.gb.gz` | present (24 files) | RefSeq records of chr1–22, X (`NC_000023`), Y (`NC_000024`), assembly GRCh38.p12 |

Note that `GB/` has **no chrM record** (`NC_012920`): mitochondrial sites get no annotation.

The three remaining files are built from those, or downloaded:

```bash
REF=/home/h/OneDrive/Unil/data/b38

# 1. b38.fa -- uncompressed genome, needed by pileup2var and samtools stats (~3 GB)
gzip -dc $REF/b38.fa.gz > $REF/b38.fa
samtools faidx $REF/b38.fa
samtools faidx $REF/b38.fa.gz        # .fai + .gzi, needed by samtools mpileup on the bgzip file

# 2. hs_rrna.fa -- human rRNA sequences (45S/18S/5.8S/28S precursor + 5S)
#    NCBI Nucleotide: NR_046235 (45S precursor, RNA45SN1) and NR_023363 (5S, RNA5S1)
efetch -db nuccore -id NR_046235.3,NR_023363.1 -format fasta > $REF/hs_rrna.fa   # EDirect
#    or, without EDirect, the same two accessions from https://www.ncbi.nlm.nih.gov/nuccore/

# 3. hisat-3n indexes (A,G three-letter mode); the genome one takes hours and a lot of RAM
#    (the repeat index is what --repeat --repeat-limit 1000 needs at mapping time)
mkdir -p $REF/HS_rep
~/src/github/hisat-3n/hisat-3n-build --base-change A,G -p 24 $REF/hs_rrna.fa $REF/HS_rep/hs_rrna
~/src/github/hisat-3n/hisat-3n-build --base-change A,G -p 24 --repeat-index $REF/b38.fa $REF/HS_rep/b38
```

Where the reference files come from, if the copy has to be rebuilt from scratch:

| What | Source |
|------|--------|
| Genome FASTA | UCSC — `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz` (matches the `chr*` naming; bgzip it: `gzip -dc hg38.fa.gz \| bgzip -@8 > b38.fa.gz`) |
| GenBank records | Per chromosome, as in `GB/`: `efetch -db nuccore -id NC_000001 -format gbwithparts \| gzip > GB/NC_000001.gb.gz` (repeat for `NC_000002` … `NC_000024`). The whole assembly in one file is also on the NCBI FTP: `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.gbff.gz` |
| rRNA sequences | NCBI Nucleotide, accessions `NR_046235` (45S precursor) and `NR_023363` (5S) |
| hisat-3n | `https://github.com/DaehwanKimLab/hisat2/tree/hisat-3n` (build with `make hisat-3n`) |
| `pileup2var` | `https://github.com/shunliubio/pileup2var` |

The chromosome names must stay consistent across the three sources: `annotPileupRev` builds the
region name as `"chr" . <chromosome of the DEFINITION line>` and uses it both for
`samtools mpileup -r` and as the FASTA identifier. The current pair (UCSC `chr1` FASTA + RefSeq
`NC_000001` records whose DEFINITION reads `Homo sapiens chromosome 1, GRCh38.p12 …`) matches;
an Ensembl FASTA (`>1`) would not.

For another organism, point `REF_DIR` at a folder with the same five entries and set
`GENOME_INDEX_NAME` accordingly (it is used in every output file name, e.g. `dm6`).

### References/Accessory files

Annotation comes from the RefSeq GenBank records, not from a GTF. Only the records whose
`DEFINITION` line matches `… chromosome <name>[.,]` are processed: **unplaced scaffolds and
records without that wording are silently skipped**. The features used are `CDS`, `mRNA`,
`misc_RNA`, `ncRNA`, `precursor_RNA`, `rRNA` and `tRNA`; `CDS` takes precedence over the others
for a position annotated twice, and features described with `order(...)`, or with a
`complement(...)` inside a `join(...)`, are skipped.

## Usage/Step-wise Execution

Everything is launched from the dataset folder, which holds `run`,
`2_map_and_count_stitched_reads.sh` and the input `fq/`, with `annotPileupRev` one level above
(see *Expected layout at run time*):

```bash
mkdir -p /path/to/work/dataset
cp /home/h/PycharmProjects/Unil/PreProcess/annotPileupRev /path/to/work/
cp /home/h/PycharmProjects/Unil/PreProcess/{run,2_map_and_count_stitched_reads.sh} /path/to/work/dataset/
cd /path/to/work/dataset      # then put the FASTQ files in fq/
chmod +x run 2_map_and_count_stitched_reads.sh ../annotPileupRev

./run                                            # default REF_DIR
./run /home/h/OneDrive/Unil/data/b38             # explicit reference folder
REF_DIR=/data6/b38 SAMPLES="S2R Y14a" ./run      # other reference folder and sample list
```

The first line of the log recalls which reference folder is in use; if a file is missing there,
both scripts stop immediately and name it.

**Step 1: Adapter trimming, round 1** — `cutadapt` removes the 5′ adapter of R1 (`-g`) and the
3′ adapter of R2 (`-A`), quality-trims at Q6 and requires R2 ≥ 46 nt. All the lanes of the
sample are concatenated on the fly by process substitution.

**Step 2: Adapter trimming, round 2** — removes the 3′ adapter of R1 (UMI + adapter, `-O 7`)
and the 5′ adapter of R2.

**Step 3: UMI extraction** — `umi_tools extract` moves the first 6 nt of R2 into the read name
of both mates.

**Step 4: Stitching** — `NGmerge` merges the overlapping mates into a single consensus read
`<sample>.stitched.fq.gz`; the non-stitched pairs are written to `<sample>.nonstitched_?.fastq.gz`
(and deleted at the end of the iteration). The trimmed intermediate FASTQ files are removed.

**Step 5: rRNA pass** — `2_map_…sh` maps the stitched reads on the rRNA index with
`hisat-3n --base-change A,G --norc --no-spliced-alignment --no-softclip`; unmapped reads
(`--un-gz`) continue to the genome, and the rRNA alignments are kept for QC (stats, conversion
filter, `pileup2var`).

**Step 6: Genome pass** — `hisat-3n --base-change A,G --repeat --repeat-limit 1000
--bowtie2-dp 0 --rna-strandness R` on the genome index, then MAPQ ≥ 60 to drop the
multi-mappers, sorting and indexing.

**Step 7: Deduplication** — `umi_tools dedup --method=unique --spliced-is-unique` collapses the
PCR duplicates; `samtools stats` is run on the deduplicated BAM.

**Step 8: Conversion filter** — reads whose converted-A fraction `Yf/(Yf+Zf)` is below
`A2G_percent` are discarded, giving the final
`outstitch/<sample>/<sample>.b38.align.sorted.dedup.flt.bam`. `pileup2var` (flag filter `-f 524`,
`-a A`, `-s R`, `-c 1`) writes the per-position conversion table. The intermediate BAMs are
deleted, and the total runtime is printed.

**Step 9: Annotation** — `gzip -dc $GB_DIR/*.gb.gz | $ANNOT_PL <bam> $GENOME_FA_GZ`.
For each chromosome the script loads the sequence, walks every transcript to tag each A with its
feature kind, amino acid, frame, 7 nt context, distance to the donor/acceptor site and gene, then
runs two `samtools mpileup` passes (common filter `--ff UNMAP,SECONDARY,QCFAIL,DUP`):

| Pass | mpileup option | Reference base kept | Counted as `nbA` / `nbG` | Annotates |
|------|----------------|---------------------|--------------------------|-----------|
| Genes on the **+** strand | `--rf REVERSE` (reverse reads only) | `A` | `,` / `g` | `+` strand table |
| Genes on the **−** strand | `--ff …,REVERSE` (forward reads only) | `T` | `.` / `C` | `−` strand table |

This split is what makes the counting strand-aware on a reverse-stranded library; a read base
appearing on the unexpected strand makes the script `die`, which is an intentional sanity check.
Sites with `nbA + nbG < 10` are not printed.

## Outputs

| Output File | Format | Script | Description |
|-------------|--------|--------|-------------|
| `<sample>.annot.txt` | TSV with header | `annotPileupRev` | **Main output**: one line per covered A site, with counts and annotation |
| `outstitch/<sample>/<sample>.b38.align.sorted.dedup.flt.bam` | BAM | `2_map_…sh` | Final alignments (unique, deduplicated, conversion-filtered) |
| `outstitch/<sample>/<sample>.b38.pileup2var.flt.txt` | TSV | `2_map_…sh` | Per-position A/G counts on the genome |
| `outstitch/<sample>/<sample>.hs_rrna.align.sorted.flt.bam` | BAM | `2_map_…sh` | rRNA alignments after conversion filter |
| `outstitch/<sample>/<sample>.hs_rrna.pileup2var.flt.txt` | TSV | `2_map_…sh` | Per-position A/G counts on the rRNA |
| `outstitch/<sample>/*.stats` | samtools stats | `2_map_…sh` | QC of the rRNA and genome alignments |
| `outstitch/<sample>/*.dedup.log`, `*_per_umi.tsv`, `*_edit_distance.tsv` | Text / TSV | `2_map_…sh` | `umi_tools dedup` log and statistics |
| `<sample>.a3a5trim.umi.round1.log` / `.round2.log` | Text | `run` | cutadapt reports |
| `<sample>.a3a5trim.umi.log` | Text | `run` | `umi_tools extract` report |
| `<sample>.stitched.fq.gz` | gzipped FASTQ | `run` | Stitched reads given to the mapper |

Columns of `<sample>.annot.txt`:

| Column | Description |
|--------|-------------|
| `chr`, `pos` | Chromosome and 1-based position |
| `nt` | Reference base: `A` = site on the **+** strand, `T` = site on the **−** strand |
| `nbA` | Reads still reading A (protected → methylated) |
| `nbG` | Reads converted to G (unmethylated) |
| `nbA+G` | Total coverage at the site (≥ 10 by construction) |
| `nbStart`, `nbStop` | Number of reads starting (`^`) or ending (`$`) at the position |
| `kind` | Feature: `CDS`, `mRNA`, `ncRNA`, `misc_RNA`, `precursor_RNA`, `rRNA`, `tRNA`, or `NA` |
| `aa`, `frame` | For a `CDS`: amino acid of the codon and position of the A in it (0, 1 or 2) |
| `context` | 7 nt window centred on the A, **transcript-oriented** (reverse-complemented on the − strand) |
| `donor`, `acceptor` | Distance (< 3 nt) to the downstream / upstream exon boundary, else `NA` |
| `gene` | Gene symbol from the GenBank `/gene=` qualifier |

These are exactly the columns expected by the downstream R pipelines (`chr`, `pos`, `nt`,
`context`, `nbA`, `nbG`, `nbA.G`, `kind`, `gene`).

