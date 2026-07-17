version 1.0

workflow NativeMismatchRate {

    input {
        Array[Array[File]] pysamstats_tsvs_per_sample
        Array[String]      sample_labels

        Array[String] chroms = [
            "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
            "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
            "chr17","chr18","chr19","chr20","chr21","chr22"
        ]

        String chroms_str = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

        File  analysis_regions_bed
        File  variants_bed
        File  junctions_bed
        File  transcript_strand_bed
        File  giab_bed
        File? ivt_summary
        File? modkit_bed

        Int   min_coverage     = 10
        Int   junction_bp      = 50
        Float min_mod_fraction = 0.5
        Int   min_mod_coverage = 20

        Int cpus       = 4
        Int memory_gb  = 40
    }

    Array[Pair[String, Array[File]]] samples = zip(sample_labels, pysamstats_tsvs_per_sample)

    scatter (sample in samples) {
        String      label      = sample.left
        Array[File] chrom_tsvs = sample.right

        Array[Pair[String, File]] chrom_pairs = zip(chroms, chrom_tsvs)

        scatter (chrom_pair in chrom_pairs) {
            String chrom     = chrom_pair.left
            File   chrom_tsv = chrom_pair.right

            call ChromTask {
                input:
                    pysamstats_tsv       = chrom_tsv,
                    chrom                = chrom,
                    sample_label         = label,
                    analysis_regions_bed = analysis_regions_bed,
                    variants_bed         = variants_bed,
                    junctions_bed        = junctions_bed,
                    transcript_strand_bed = transcript_strand_bed,
                    min_coverage         = min_coverage,
                    junction_bp          = junction_bp,
                    cpus                 = cpus,
                    memory_gb            = memory_gb
            }
        }

        # Collect per-chrom position TSVs for merging
        Array[File] valid_position_tsvs = ChromTask.per_position_tsv

        call MergeTask {
            input:
                position_tsvs        = valid_position_tsvs,
                sample_label         = label,
                analysis_regions_bed = analysis_regions_bed,
                variants_bed         = variants_bed,
                junctions_bed        = junctions_bed,
                giab_bed             = giab_bed,
                ivt_summary          = ivt_summary,
                modkit_bed           = modkit_bed,
                min_coverage         = min_coverage,
                junction_bp          = junction_bp,
                min_mod_fraction     = min_mod_fraction,
                min_mod_coverage     = min_mod_coverage,
                chroms_str           = chroms_str,
                cpus                 = cpus,
                memory_gb            = memory_gb
        }
    }

    output {
        Array[File]        mismatch_summaries   = MergeTask.mismatch_summary
        Array[File]        masking_reports      = MergeTask.masking_report
        Array[File]        per_position_tsvs    = MergeTask.per_position_tsv
        Array[File?]       substitution_spectra = MergeTask.substitution_spectrum
        Array[File?]       modkit_stratified    = MergeTask.modkit_stratified
        Array[File?]       native_ivt_deltas    = MergeTask.native_ivt_delta
        Array[Array[File]] chrom_position_tsvs  = ChromTask.per_position_tsv
    }
}


task ChromTask {

    input {
        File   pysamstats_tsv
        String chrom
        String sample_label

        File   analysis_regions_bed
        File   variants_bed
        File   junctions_bed
        File   transcript_strand_bed

        Int    min_coverage
        Int    junction_bp

        Int    cpus
        Int    memory_gb
    }

    command <<<
        set -euo pipefail

        echo "[WDL ChromTask] ~{sample_label} / ~{chrom}"

        mkdir -p masks_dir outdir

        ln -sf ~{analysis_regions_bed}   masks_dir/analysis_regions.bed
        ln -sf ~{variants_bed}           masks_dir/variants.bed
        ln -sf ~{junctions_bed}          masks_dir/junctions.bed
        ln -sf ~{transcript_strand_bed}  masks_dir/transcript_strand.bed

        python /usr/local/bin/calculate_native_mismatch_rates_fast.py \
            --mode         chrom             \
            --pysamstats   ~{pysamstats_tsv} \
            --masks-dir    masks_dir/        \
            --chrom        ~{chrom}          \
            --sample-label ~{sample_label}   \
            --min-coverage ~{min_coverage}   \
            --junction-bp  ~{junction_bp}    \
            --outdir       outdir/

        echo "[WDL ChromTask] Done: ~{sample_label} / ~{chrom}"
    >>>

    output {
        File per_position_tsv = "outdir/~{sample_label}_~{chrom}_per_position_mismatch.tsv"
    }

    runtime {
        cpu:             cpus
        memory:          "~{memory_gb}GB"
        maxRunTime:      21600
        runtime_minutes: 240
        docker:          "kandarprj/drs-native-mismatch-1.5"
    }
}


task MergeTask {

    input {
        Array[File] position_tsvs
        String      sample_label

        File   analysis_regions_bed
        File   variants_bed
        File   junctions_bed
        File   giab_bed
        File?  ivt_summary
        File?  modkit_bed

        Int    min_coverage
        Int    junction_bp
        Float  min_mod_fraction
        Int    min_mod_coverage
        String chroms_str

        Int    cpus
        Int    memory_gb
    }

    command <<<
        set -euo pipefail

        echo "[WDL MergeTask] ~{sample_label} — ~{length(position_tsvs)} chrom TSVs"

        mkdir -p masks_dir outdir

        ln -sf ~{analysis_regions_bed} masks_dir/analysis_regions.bed
        ln -sf ~{variants_bed}         masks_dir/variants.bed
        ln -sf ~{junctions_bed}        masks_dir/junctions.bed

        VDIR=$(dirname ~{variants_bed})
        [ -f "$VDIR/giab_minus_variants.bed" ] && \
            ln -sf "$VDIR/giab_minus_variants.bed" masks_dir/giab_minus_variants.bed || true
        [ -f "$VDIR/minus_junctions.bed" ] && \
            ln -sf "$VDIR/minus_junctions.bed" masks_dir/minus_junctions.bed || true

        MODKIT_ARG=""
        if [ -n "~{default="" modkit_bed}" ]; then
            MODKIT_ARG="--modkit ~{modkit_bed}"
        fi

        IVT_ARG=""
        if [ -n "~{default="" ivt_summary}" ]; then
            IVT_ARG="--ivt-summary ~{ivt_summary}"
        fi

        python /usr/local/bin/calculate_native_mismatch_rates_fast.py \
            --mode             merge                     \
            --position-tsvs    ~{sep=" " position_tsvs} \
            --masks-dir-merge  masks_dir/               \
            --giab-bed         ~{giab_bed}               \
            --sample-label     ~{sample_label}           \
            --min-coverage     ~{min_coverage}           \
            --junction-bp      ~{junction_bp}            \
            --chroms           ~{chroms_str}             \
            --outdir           outdir/                   \
            ${MODKIT_ARG}                                \
            ${IVT_ARG}

        echo "[WDL MergeTask] Done: ~{sample_label}"
    >>>

    output {
        File  mismatch_summary      = "outdir/~{sample_label}_mismatch_rate_summary.tsv"
        File  masking_report        = "outdir/~{sample_label}_masking_report.tsv"
        File  per_position_tsv      = "outdir/~{sample_label}_per_position_mismatch.tsv"
        File? substitution_spectrum = "outdir/~{sample_label}_substitution_spectrum.tsv"
        File? modkit_stratified     = "outdir/~{sample_label}_modkit_stratified_rates.tsv"
        File? native_ivt_delta      = "outdir/~{sample_label}_native_ivt_delta.tsv"
    }

    runtime {
        cpu:             cpus
        memory:          "~{memory_gb}GB"
        maxRunTime:      21600
        runtime_minutes: 240
        docker:          "kandarprj/drs-native-mismatch-1.5"
    }
}
