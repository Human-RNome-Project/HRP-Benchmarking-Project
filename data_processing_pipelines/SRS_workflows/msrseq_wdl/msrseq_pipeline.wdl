version 1.0

struct MSRSample {
    String sample_id
    String treatment
    Int replicate
    String barcode
    File fastq_read2
}

workflow msrseq_pipeline {
    input {
        File mature_trna_fasta
        File chromosomal_trna_fasta
        File isodecoder_table
        Array[MSRSample] samples
        
        String docker_image = "msrseq-pipeline:v1"
        
        Int ref_cpus = 8
        String ref_memory = "16GB"
        Int ref_runtime_minutes = 30
        
        Int sample_cpus = 16
        String sample_memory = "64GB"
        Int sample_runtime_minutes = 120
        
        Int agg_cpus = 4
        String agg_memory = "16GB"
        Int agg_runtime_minutes = 30
        
        Int convert_cpus = 4
        String convert_memory = "16GB"
        Int convert_runtime_minutes = 30
    }

    call PrepareReference {
        input:
            mature_trna_fasta = mature_trna_fasta,
            docker_image = docker_image,
            cpus = ref_cpus,
            memory = ref_memory,
            runtime_minutes = ref_runtime_minutes
    }

    scatter (s in samples) {
        call ProcessSample {
            input:
                sample = s,
                cleaned_fasta = PrepareReference.cleaned_fasta,
                bt2_index_files = PrepareReference.bt2_index_files,
                bt2_index_prefix = PrepareReference.bt2_index_prefix,
                docker_image = docker_image,
                cpus = sample_cpus,
                memory = sample_memory,
                runtime_minutes = sample_runtime_minutes
        }
    }

    call AggregateTsvs {
        input:
            all_tsvs = flatten(ProcessSample.sample_tsvs),
            docker_image = docker_image,
            cpus = agg_cpus,
            memory = agg_memory,
            runtime_minutes = agg_runtime_minutes
    }

    call BedrmodConversion {
        input:
            data_cleaned_zip = AggregateTsvs.data_cleaned_zip,
            chromosomal_trna_fasta = chromosomal_trna_fasta,
            isodecoder_table = isodecoder_table,
            docker_image = docker_image,
            cpus = convert_cpus,
            memory = convert_memory,
            runtime_minutes = convert_runtime_minutes
    }

    output {
        File cleaned_ref_fasta        = PrepareReference.cleaned_fasta
        Array[File] bt2_index_files   = PrepareReference.bt2_index_files
        Array[File] all_tsvs          = flatten(ProcessSample.sample_tsvs)
        File data_cleaned_zip         = AggregateTsvs.data_cleaned_zip
        Array[File] bedrmod_files     = BedrmodConversion.bedrmod_files
        File bs_rep1                  = BedrmodConversion.bs_rep1
        File bs_rep2                  = BedrmodConversion.bs_rep2
        File bs_rep3                  = BedrmodConversion.bs_rep3
        File cbh_rep1                 = BedrmodConversion.cbh_rep1
        File cbh_rep2                 = BedrmodConversion.cbh_rep2
        File cbh_rep3                 = BedrmodConversion.cbh_rep3
        File ctrl_mut_rep1            = BedrmodConversion.ctrl_mut_rep1
        File ctrl_mut_rep2            = BedrmodConversion.ctrl_mut_rep2
        File ctrl_mut_rep3            = BedrmodConversion.ctrl_mut_rep3
        File ctrl_del_rep1            = BedrmodConversion.ctrl_del_rep1
        File ctrl_del_rep2            = BedrmodConversion.ctrl_del_rep2
        File ctrl_del_rep3            = BedrmodConversion.ctrl_del_rep3
    }
}

task PrepareReference {
    input {
        File mature_trna_fasta
        String docker_image
        Int cpus
        String memory
        Int runtime_minutes
    }

    command <<<
        set -euo pipefail
        
        Rscript /usr/local/bin/refseq_cleanup.R "~{mature_trna_fasta}" "mature_tRNAs_cleaned.fa"
        
        mkdir -p bt2_index
        bowtie2-build --threads ~{cpus} -f "mature_tRNAs_cleaned.fa" "bt2_index/mature_tRNAs_cleaned"
    >>>

    output {
        File cleaned_fasta         = "mature_tRNAs_cleaned.fa"
        Array[File] bt2_index_files = glob("bt2_index/*")
        String bt2_index_prefix    = "bt2_index/mature_tRNAs_cleaned"
    }

    runtime {
        cpu: cpus
        memory: memory
        maxRunTime: runtime_minutes * 60
        runtime_minutes: runtime_minutes
        docker: docker_image
    }
}

task ProcessSample {
    input {
        MSRSample sample
        File cleaned_fasta
        Array[File] bt2_index_files
        String bt2_index_prefix
        String docker_image
        Int cpus
        String memory
        Int runtime_minutes
    }

    command <<<
        set -euo pipefail
        
        mkdir -p sam_out bam_out tsv_out
        
        # Link bowtie2 index files to preserve index prefix path
        idx_dir=$(dirname "~{bt2_index_prefix}")
        mkdir -p "${idx_dir}"
        for f in ~{sep=' ' bt2_index_files}; do
            ln -sf "$f" "${idx_dir}/$(basename "$f")"
        done
        
        # 1. Bowtie2 sense Read 2 local alignment
        bowtie2 -x "~{bt2_index_prefix}" -U "~{sample.fastq_read2}" -S "sam_out/~{sample.sample_id}.sam" -q -p ~{cpus} --local --no-unal
        
        # 2. Extension-length SAM binning
        python3 /usr/local/bin/sam_bin_split.py -i "sam_out/~{sample.sample_id}.sam" -o "sam_out/" -breaks 0,10,20,30,40,50,60
        
        # 3. For each non-empty SAM: BAM conversion, sorting, IGV coverage count, and base-level TSV conversion
        for sam in sam_out/~{sample.sample_id}*.sam; do
            # Skip empty header-only SAM files with 0 alignments
            if ! grep -q -m 1 "^[^@]" "$sam"; then
                continue
            fi
            
            base=$(basename "$sam" .sam)
            bam="bam_out/${base}.bam"
            sort_bam="bam_out/${base}.sort.bam"
            wig="bam_out/${base}.wig"
            tsv="tsv_out/${base}.sam.tsv"
            
            samtools view -@ 4 -bS -o "$bam" "$sam"
            samtools sort -@ 4 "$bam" -o "$sort_bam"
            igvtools count -z 5 -w 1 -e 250 --bases "$sort_bam" "$wig" "~{cleaned_fasta}"
            python3 /usr/local/bin/wig_to_tsv_low_mem_2.py -i "$wig" -r "~{cleaned_fasta}" -o "$tsv"
            
            rm -f "$bam" "$sort_bam" "$wig"
        done
        
        rm -rf sam_out bam_out
    >>>

    output {
        Array[File] sample_tsvs = glob("tsv_out/*.tsv")
    }

    runtime {
        cpu: cpus
        memory: memory
        maxRunTime: runtime_minutes * 60
        runtime_minutes: runtime_minutes
        docker: docker_image
    }
}

task AggregateTsvs {
    input {
        Array[File] all_tsvs
        String docker_image
        Int cpus
        String memory
        Int runtime_minutes
    }

    command <<<
        set -euo pipefail
        
        mkdir -p tsv_input
        for tsv in ~{sep=' ' all_tsvs}; do
            ln -sf "$tsv" "tsv_input/$(basename "$tsv")"
        done
        
        Rscript /usr/local/bin/tRNAmod_cleaner.R "tsv_input" "data_cleaned_5_HRPC.csv" "data_cleaned_5_HRPC.csv.zip"
    >>>

    output {
        File data_cleaned_zip = "data_cleaned_5_HRPC.csv.zip"
        File data_cleaned_csv = "data_cleaned_5_HRPC.csv"
    }

    runtime {
        cpu: cpus
        memory: memory
        maxRunTime: runtime_minutes * 60
        runtime_minutes: runtime_minutes
        docker: docker_image
    }
}

task BedrmodConversion {
    input {
        File data_cleaned_zip
        File chromosomal_trna_fasta
        File isodecoder_table
        String docker_image
        Int cpus
        String memory
        Int runtime_minutes
    }

    command <<<
        set -euo pipefail
        
        mkdir -p bedrmod_out
        score_script="/opt/msrseq/scripts/bedrmod_score_stoichiometry.R"
        
        Rscript /usr/local/bin/bedrmod_converter.R \
            "~{data_cleaned_zip}" \
            "~{chromosomal_trna_fasta}" \
            "~{isodecoder_table}" \
            "${score_script}" \
            "bedrmod_out"
    >>>

    output {
        Array[File] bedrmod_files = glob("bedrmod_out/*.bedRMod")
        File bs_rep1              = "bedrmod_out/HRPC_BS_1.bedRMod"
        File bs_rep2              = "bedrmod_out/HRPC_BS_2.bedRMod"
        File bs_rep3              = "bedrmod_out/HRPC_BS_3.bedRMod"
        File cbh_rep1             = "bedrmod_out/HRPC_CBH_1.bedRMod"
        File cbh_rep2             = "bedrmod_out/HRPC_CBH_2.bedRMod"
        File cbh_rep3             = "bedrmod_out/HRPC_CBH_3.bedRMod"
        File ctrl_mut_rep1        = "bedrmod_out/HRPC_ctrl_mut_1.bedRMod"
        File ctrl_mut_rep2        = "bedrmod_out/HRPC_ctrl_mut_2.bedRMod"
        File ctrl_mut_rep3        = "bedrmod_out/HRPC_ctrl_mut_3.bedRMod"
        File ctrl_del_rep1        = "bedrmod_out/HRPC_ctrl_del_1.bedRMod"
        File ctrl_del_rep2        = "bedrmod_out/HRPC_ctrl_del_2.bedRMod"
        File ctrl_del_rep3        = "bedrmod_out/HRPC_ctrl_del_3.bedRMod"
    }

    runtime {
        cpu: cpus
        memory: memory
        maxRunTime: runtime_minutes * 60
        runtime_minutes: runtime_minutes
        docker: docker_image
    }
}
