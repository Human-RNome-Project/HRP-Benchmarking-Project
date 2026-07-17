version 1.0

## =============================================================================
## pipeline_read_error_rates.wdl
##
## Per-read DRS error rate calculation — parallel by chromosome.
##
## Scatter: one ChromTask per chromosome, each sampling --sample-n reads
##          from that chrom of the BAM using reservoir sampling.
## Merge:   MergeTask concatenates per-read TSVs, aggregates count arrays,
##          and produces final summary outputs.
##
## Total reads sampled = n_chroms × sample_n_per_chrom
## (e.g. 22 × 2274 ≈ 50,000 reads genome-wide)
## =============================================================================

workflow ReadErrorRates {

    input {
        File          bam
        File          bai
        String        sample_label

        Array[String] chroms = [
            "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
            "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
            "chr17","chr18","chr19","chr20","chr21","chr22",
            "chrX","chrY","chrM"
        ]

        # Reads sampled per chromosome
        # Set to ceil(total_sample_n / n_chroms): e.g. 50000/25 = 2000
        Int   sample_n_per_chrom
        Int   min_mapq
        Int   min_length
        Int   seed

        Int   cpus
    }

    scatter (chrom in chroms) {
        call ChromTask {
            input:
                bam          = bam,
                bai          = bai,
                chrom        = chrom,
                sample_label = sample_label,
                sample_n     = sample_n_per_chrom,
                min_mapq     = min_mapq,
                min_length   = min_length,
                seed         = seed,
                cpus         = cpus
        }
    }

    call MergeTask {
        input:
            per_read_tsvs = ChromTask.per_read_tsv,
            npz_files     = ChromTask.counts_npz,
            sample_label  = sample_label,
            sample_n      = sample_n_per_chrom * length(chroms),
            seed          = seed,
            cpus          = cpus
    }

    output {
        Array[File] chrom_per_read_tsvs = ChromTask.per_read_tsv
        Array[File] chrom_counts_npzs   = ChromTask.counts_npz
        File        per_read_tsv        = MergeTask.per_read_tsv
        File        summary_tsv         = MergeTask.summary_tsv
        File        spectrum_tsv        = MergeTask.spectrum_tsv
        File        positional_bias_tsv = MergeTask.positional_bias_tsv
    }

    meta {
        author      : "K"
        description : "Per-read DRS error rates — parallel chromosome scatter"
    }
}


task ChromTask {

    input {
        File   bam
        File   bai
        String chrom
        String sample_label
        Int    sample_n
        Int    min_mapq
        Int    min_length
        Int    seed
        Int    cpus
    }

    command <<<
        set -euo pipefail

        echo "[WDL ChromTask] ~{sample_label} / ~{chrom}"

        # BAI must be co-located with BAM for pysam
        BAM_LOCAL=$(basename ~{bam})
        BAI_LOCAL="${BAM_LOCAL}.bai"
        ln -sf ~{bam} "$BAM_LOCAL"
        ln -sf ~{bai} "$BAI_LOCAL"

        mkdir -p outdir

        python /usr/local/bin/calculate_read_error_rates_parallel.py \
            --mode         chrom              \
            --bam          "$BAM_LOCAL"       \
            --chrom        ~{chrom}           \
            --sample-n     ~{sample_n}        \
            --min-mapq     ~{min_mapq}        \
            --min-length   ~{min_length}      \
            --sample-label ~{sample_label}    \
            --seed         ~{seed}            \
            --outdir       outdir/

        echo "[WDL ChromTask] Done: ~{sample_label} / ~{chrom}"
    >>>

    output {
        File per_read_tsv = "outdir/~{sample_label}_~{chrom}_per_read.tsv"
        File counts_npz   = "outdir/~{sample_label}_~{chrom}_counts.npz"
    }

    runtime {
        cpu:             cpus
        memory:          "40GB"
        maxRunTime:      7200
        runtime_minutes: 120
        docker:          "kandarprj/drs-read-error-rates-1.0"
    }
}


task MergeTask {

    input {
        Array[File] per_read_tsvs
        Array[File] npz_files
        String      sample_label
        Int         sample_n
        Int         seed
        Int         cpus
    }

    command <<<
        set -euo pipefail

        echo "[WDL MergeTask] ~{sample_label} — merging ~{length(per_read_tsvs)} chroms"

        mkdir -p outdir

        python /usr/local/bin/calculate_read_error_rates_parallel.py \
            --mode           merge                        \
            --per-read-tsvs  ~{sep=" " per_read_tsvs}   \
            --npz-files      ~{sep=" " npz_files}        \
            --sample-n       ~{sample_n}                 \
            --sample-label   ~{sample_label}             \
            --seed           ~{seed}                     \
            --outdir         outdir/

        echo "[WDL MergeTask] Done: ~{sample_label}"
    >>>

    output {
        File per_read_tsv        = "outdir/~{sample_label}_per_read_error_rates.tsv"
        File summary_tsv         = "outdir/~{sample_label}_read_error_rate_summary.tsv"
        File spectrum_tsv        = "outdir/~{sample_label}_read_substitution_spectrum.tsv"
        File positional_bias_tsv = "outdir/~{sample_label}_read_positional_bias.tsv"
    }

    runtime {
        cpu:             cpus
        memory:          "40GB"
        maxRunTime:      7200
        runtime_minutes: 120
        docker:          "kandarprj/drs-read-error-rates-1.0"
    }
}
