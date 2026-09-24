# Processing of UMBS-seq (UBS-seq) Short-Read Sequencing Data (m5C_UMBS_seq_rep1 / rep2)

End-to-end workflow for the m5C RNome samples of the Human RNome Project, from raw FASTQ
files to a bedRMod file.

---

# Data Download Instructions

1. Open the project's data repository and follow the instructions under "Raw data".
2. Download the paired-end Illumina NovaSeq X FASTQ files from the "Short-read Sequencing" folder below into `${rawdir}`.

| HRP sample name         | Internal ID | Raw FASTQ files                                                        |
| ----------------------- | ----------- | ---------------------------------------------------------------------- |
| `m5C_UMBS_seq_rep1`  | `YSL-5`     | `HRP_B_032_1_R1.fastq.gz`, `HRP_B_032_1_R2.fastq.gz` |
| `m5C_UMBS_seq_rep2`  | `YSL-6`     | `HRP_B_032_2_R1.fastq.gz`, `HRP_B_032_2_R2.fastq.gz` |

The two libraries are technical replicates of UMBS-treated (bisulfite-converted) human mRNA.
Throughout the pipeline the samples are referred to by their internal IDs `YSL-5` and `YSL-6`,
which is also how they are named in the intermediate files.

---

# Requirements & Setup

A conda environment file (`UMBS_pipeline.yml`) is provided.

```sh
conda env create -f UMBS_pipeline.yml
conda activate UMBS_pipeline
```

| Software        | Version used | Purpose                                     |
| --------------- | ------------ | ------------------------------------------- |
| cutadapt        | 5.2          | adapter / quality trimming                  |
| hisat-3n        | 0.0.3 (hisat2 2.2.1) | C-to-T aware alignment              |
| hisat-3n-table  | 0.0.3        | per-position conversion counts              |
| samtools        | 1.23.1       | BAM handling, tag-based read filtering      |
| htslib / bgzip  | 1.23.1       | compressed TSV output                       |
| GATK            | 4.5          | duplicate marking                           |
| python          | >= 3.9       | motif annotation, m5C calling, bedRMod      |
| numpy / pandas / scipy | any recent | m5C calling                            |

`hisat-3n` and `hisat-3n-table` come from the `hisat-3n` bioconda package; the plain
`hisat2` package does not contain them.

---

# Input

- Paired-end FASTQ files, no pre-processing required.
- Read 1 and read 2 both start with 5 nt that are removed during trimming (`-u 5 -U 5`);
  the removed prefixes are appended to the read name so that they remain available.
- The library is directional C-to-T (`--base-change C,T`).

---

# Workflow

## 1. Prepare the reference and the hisat-3n index

The alignment reference is the human genome combined with the spike-in controls that are
added to every UBS-seq library (unmethylated lambda phage DNA-derived RNA and pUC19-derived
RNA, used to measure the conversion rate). Reads are mapped once against the combined
reference and split by reference afterwards.

```sh
refdir='your_reference_directory'

# human genome (Ensembl GRCh38 primary assembly)
curl -O https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > ${refdir}/human.fa

# combined reference (human + spike-ins)
cat ${refdir}/human.fa ${refdir}/lambda.fa ${refdir}/pUC19.fa > ${refdir}/refs_gen.fa
samtools faidx ${refdir}/refs_gen.fa

# C-to-T index
hisat-3n-build --base-change C,T -p 32 ${refdir}/refs_gen.fa ${refdir}/hisat3n/refs_gen
```

If the spike-ins are not needed, the same commands with `human.fa` alone reproduce the human
calls; only the conversion-rate QC is then unavailable.

## 2. Read processing, alignment and per-position conversion counts

All steps are implemented in `pipeline.sh` (provided next to this README); it reproduces the
core pipeline of the original analysis notebook `02_core_pipeline.ipynb`.

```sh
./pipeline.sh --index ${refdir}/hisat3n/refs_gen \
              --ref   ${refdir}/refs_gen.fa \
              --threads 32
```

The individual steps and their parameters:

### trim (cutadapt 5.2)

```sh
cutadapt -j 0 -n 2 \
    -a 'AGATCGGAAGAGCACACGTCT;e=0.15;o=4;anywhere;' \
    -A 'AGATCGGAAGAGCGTCGTGT;e=0.15;o=4;anywhere;' \
    --max-n=0 -q 15 --nextseq-trim=15 -m 20 \
    -u 5 -U 5 --rename='{id}_{r1.cut_prefix}{r2.cut_prefix}' \
    --too-short-output=trim/${SAMPLE}_R1.fa_short \
    --too-short-paired-output=trim/${SAMPLE}_R2.fa_short \
    -o trim/${SAMPLE}_R1.fq.gz -p trim/${SAMPLE}_R2.fq.gz \
    ${SAMPLE}_R1.fq.gz ${SAMPLE}_R2.fq.gz > trim/${SAMPLE}.report
```

### map (hisat-3n 0.0.3)

```sh
hisat-3n --index ${INDEX} -p ${THREADS} --base-change C,T --mp 8,2 \
    --no-spliced-alignment \
    --summary-file map/${SAMPLE}.summary --new-summary \
    -1 trim/${SAMPLE}_R1.fq.gz -2 trim/${SAMPLE}_R2.fq.gz \
    --un-conc-gz map/${SAMPLE}_R%.fq.gz -S map/${SAMPLE}.sam
samtools view -@ ${THREADS} -F4 -b map/${SAMPLE}.sam \
  | samtools sort -@ ${THREADS} --write-index -O BAM -o map/${SAMPLE}.bam -
```

### mark_duplicates and dedup (GATK 4.5 + samtools)

```sh
gatk AddOrReplaceReadGroups -I map/${SAMPLE}.bam -O mark_duplicates/${SAMPLE}.rg.bam \
    -LB 1 -PL ILLUMINA -PU 1 -SM ${SAMPLE}
gatk MarkDuplicates -I mark_duplicates/${SAMPLE}.rg.bam -O mark_duplicates/${SAMPLE}.bam \
    -M mark_duplicates/${SAMPLE}.tmp \
    --ASSUME_SORT_ORDER coordinate --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500
gatk CollectDuplicateMetrics -I mark_duplicates/${SAMPLE}.bam -M mark_duplicates/${SAMPLE}.metrics
samtools view -h -F 1024 mark_duplicates/${SAMPLE}.bam -o dedup/${SAMPLE}.bam
```

### split by reference

```sh
samtools view -b dedup/${SAMPLE}.bam ${HUMAN_CONTIGS} -o split/${SAMPLE}.human.bam
```

(the same is done for `lambda` and `pUC19`, which are used only for conversion-rate QC)

### dedup_filter - read-level quality filter

Reads are kept only when all three conditions hold, using the hisat-3n tags
`Yf` (converted Cs), `Zf` (unconverted Cs) and `XM` (mismatches):

```sh
samtools view -e '[XM] * 20 <= (qlen - sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]' \
    -O BAM -o dedup_filter/${SAMPLE}.bam split/${SAMPLE}.human.bam
```

- `[XM] * 20 <= (qlen - sclen)` - fewer than one mismatch per 20 aligned bases
- `[Zf] <= 3` - at most 3 unconverted Cs per read
- `3 * [Zf] <= [Zf] + [Yf]` - unconverted Cs are less than one third of all Cs of the read

### conv_unconv3n_filter - per-position table

```sh
samtools view -e "rlen<100000" -h dedup_filter/${SAMPLE}.bam \
  | hisat-3n-table -p ${THREADS} --unique-only --alignments - \
        --ref ${REF_FASTA} --output-name /dev/stdout --base-change C,T \
  | python annotate_motif.py --ref ${REF_FASTA} --sample ${SAMPLE} \
  | bgzip -@ ${THREADS} -c > conv_unconv3n_filter/${SAMPLE}.tsv.gz
```

`annotate_motif.py` reformats the raw `hisat-3n-table` output and appends the trinucleotide
context, producing the table that the m5C calling consumes:

| Column | Meaning |
| ------ | ------- |
| `Sample` | sample ID |
| `Chrom`, `Pos`, `Strand` | 1-based position of the cytosine |
| `Converted` | C-to-T reads |
| `Unconverted` | C-to-C reads |
| `Depth` | `Converted + Unconverted` |
| `Ratio_conv`, `Ratio_unconv` | the two ratios |
| `Motif` | trinucleotide context (CpG / CHG / CHH) |

The same table without the read-level filter (`conv_unconv3n/${SAMPLE}.tsv.gz`) is produced as
well; it is used to estimate the background non-conversion rate (next section).

## 3. m5C site calling

Implemented in `m5C_calling.py`, reproducing `05_m5C_calling.ipynb`.

At every cytosine with depth *d* and *k* unconverted reads, the probability of observing at
least *k* unconverted reads by chance is computed from the binomial tail
P(X >= *k*) ~ Binomial(*d*, *e*) (Legrand et al. 2017, as adapted by Dai et al. 2024), where
*e* is the per-sample background non-conversion rate.

**Background non-conversion rate.** *e* is estimated from the **unfiltered** per-position
table, because the read-level filter biases the filtered set towards converted reads and
would give an anti-conservative estimate:

| Sample | `e`        |
| ------ | ---------- |
| YSL-5  | `0.000731` |
| YSL-6  | `0.000836` |

**Calling thresholds**

| Parameter              | Value                                        |
| ---------------------- | -------------------------------------------- |
| minimum depth          | 3                                            |
| minimum stoichiometry  | 0 (no floor on `Ratio_unconv`)               |
| p-value threshold      | 1e-12 (stringent) **and** 1e-6 (Dai 2024 / Legrand 2017 default), run separately |

```sh
python m5C_calling.py \
    --input conv_unconv3n_filter/YSL-5.tsv.gz --sample YSL-5 \
    --error-rate 0.000731 --min-depth 3 --min-ratio 0.0 --max-p-value 1e-12 \
    --output m5C_sites_YSL-5_dpth_3_p_val_1e-12_min_rat_0.0.tsv

python m5C_calling.py \
    --input conv_unconv3n_filter/YSL-6.tsv.gz --sample YSL-6 \
    --error-rate 0.000836 --min-depth 3 --min-ratio 0.0 --max-p-value 1e-12 \
    --output m5C_sites_YSL-6_dpth_3_p_val_1e-12_min_rat_0.0.tsv

# sites called in both replicates (long format: one row per site and sample)
python m5C_calling.py --intersect \
    --input m5C_sites_YSL-5_dpth_3_p_val_1e-12_min_rat_0.0.tsv \
            m5C_sites_YSL-6_dpth_3_p_val_1e-12_min_rat_0.0.tsv \
    --output m5C_sites_intersection_YSL-5_YSL-6_dpth_3_p_val_1e-12_min_rat_0.0.tsv
```

Repeat with `--max-p-value 1e-6` for the sensitive call set.  `--error-rate` may be omitted,
in which case it is estimated from the input table itself as
`sum(Unconverted) / sum(Depth)`; passing the values above reproduces the published call sets
exactly.

The called-site tables carry 12 columns:
`Sample, Chrom, Pos, Strand, Motif, Converted, Unconverted, Depth, Ratio_conv, Ratio_unconv,
p_val, ID` (`ID` = `Chrom_Pos`, the key used for the intersection).

## 4. Convert the results into bedRMod

```sh
python convert_to_bedRmod.py \
    --input m5C_sites_intersection_YSL-5_YSL-6_dpth_3_p_val_1e-12_min_rat_0.0.tsv \
    --output m5C_UMBS_seq.bedrmod \
    --combine-samples

# per replicate
python convert_to_bedRmod.py --input m5C_sites_YSL-5_dpth_3_p_val_1e-12_min_rat_0.0.tsv \
    --output m5C_UMBS_seq_rep1.bedrmod
python convert_to_bedRmod.py --input m5C_sites_YSL-6_dpth_3_p_val_1e-12_min_rat_0.0.tsv \
    --output m5C_UMBS_seq_rep2.bedrmod
```

Conventions written by the script:

- `chromStart` is 0-based (`Pos - 1`), `chromEnd = chromStart + 1`;
- `name` is the modification short name `m5C` (`--name-field id` writes `20607`);
- `score` is `min(1000, coverage)`;
- `coverage` is `Depth` (`--combine-samples` sums the depths of the two replicates);
- `frequency` is `Ratio_unconv` as a **percentage** (0-100, two decimals);
- only called sites are written - positions that are covered but not called are not part of
  the bedRMod file.

---

# Output

- `m5C_UMBS_seq.bedrmod` - m5C sites called in both replicates (high-confidence set).
- `m5C_UMBS_seq_rep1.bedrmod`, `m5C_UMBS_seq_rep2.bedrmod` - per-replicate sites.

Header of the produced files:

```
#fileformat=bedRModv2
#organism=9606
#modification_type=RNA
#modification_names=20607:m5C:C
#assembly=GRCh38
#annotation_source=Ensembl
#annotation_version=110
#sequencing_platform=Illumina NovaSeq X
#bioinformatics_workflow=UBS-seq (hisat-3n 0.0.3) + binomial m5C calling
```

---

# References

- Dai Q. *et al.* Ultrafast bisulfite sequencing detection of 5-methylcytosine in DNA and RNA.
  *Nat. Biotechnol.* (2024).
- Legrand C. *et al.* Statistically robust methylation calling for whole-transcriptome
  bisulfite sequencing reveals distinct methylation patterns for mouse RNAs.
  *Genome Res.* 27, 1589-1596 (2017).
