version 1.0

workflow ont_tRNA_merge_modkit_grid {
    input {
        File bamfiles
        File bamindices
        String sample_id
        File reference
        Int cpus
        Array[Float] mod_thresholds
    }

    #call MergeBams {
    #    input:
    #        bams = bamfiles,
    #        bai = bamindices,
    #        sample_id = sample_id,
    #        cpus = cpus
    #}

    scatter (threshold in mod_thresholds) {
        call ModkitPileup {
            input:
                bam = bamfiles,
                bai = bamindices,
                sample_id = sample_id,
                reference = reference,
                cpus = cpus,
                threshold = threshold
        }
    }

    call NanoCompQC {
        input:
            bam = bamfiles,
            sample_id = sample_id,
            cpus = cpus
    }

    output {
#        File merged_bam = MergeBams.merged_bam
#        File merged_bai = MergeBams.merged_bai
        Array[File] modkit_beds = ModkitPileup.bed
        Array[File] modkit_logs = ModkitPileup.log
        File? nanocomp_report = NanoCompQC.report_tar
    }
}

#task MergeBams {
#    input {
#        Array[File] bams
#        Array[File] bai
#        String sample_id
#        Int cpus
#    }
#
#   command <<<
#   set -euo pipefail
#
    # Merge all aligned BAMs into a single sample BAM.
#    ls ~{sep=' ' bams} > bam_list.txt
#    samtools merge -@ ~{cpus} -o "~{sample_id}.merged.bam" -b bam_list.txt
#    samtools index -b -@ ~{cpus} "~{sample_id}.merged.bam"
#    >>>

#    output {
#        File merged_bam = "~{sample_id}.merged.bam"
#        File merged_bai = "~{sample_id}.merged.bam.bai"
#    }

#    runtime {
#        cpu: cpus
#        memory: "128GB"
#        maxRunTime: 86400 #24 hours
#        runtime_minutes: 1 #24 hours
#        docker: "nanozoo/minimap2:2.28--9e3bd01"
#    }
#}

task NanoCompQC {
    input {
        File bam
        String sample_id
        Int cpus
    }

    command <<<
    set -euo pipefail

    mkdir -p "nanocomp_report_~{sample_id}"
    filename=$(basename "~{bam}" .bam)
    NanoComp --bam "~{bam}" \
    --threads ~{cpus} \
    --outdir "nanocomp_report_${filename}"
    tar -czf "nanocomp_report_${filename}.tar.gz" "nanocomp_report_${filename}"
    >>>

    output {
        File report_tar = "nanocomp_report_" + basename(bam, ".bam") + ".tar.gz"
    }

    runtime {
        cpu: cpus
        memory: "128GB"
        maxRunTime: 43200 #12 hours (12 * 3600 seconds)
        runtime_minutes: 1 #24 hours (24 * 60 minutes)
        docker: "luxendr13/nanocomp:0.6.0"
        failOnStderr: false
    }
}

task ModkitPileup {
    input {
        File bam
        File bai
        String sample_id
        File reference
        Int cpus
        Float threshold
    }

    command <<<
    set -euo pipefail
    filename=$(basename "~{bam}" .bam)
    modkit pileup \
        ~{bam} \
        "${filename}.threshold_~{threshold}.tRNA.bed" \
        --ref ~{reference} \
        --threads ~{cpus} \
        --filter-threshold 0.8 \
        --modified-bases m5C 2OmeC inosine m6A 2OmeA pseU 2OmeU 2OmeG \
        --mod-threshold m:~{threshold} \
        --mod-threshold 19228:~{threshold} \
        --mod-threshold 17596:~{threshold} \
        --mod-threshold a:~{threshold} \
        --mod-threshold 69426:~{threshold} \
        --mod-threshold 17802:~{threshold} \
        --mod-threshold 19227:~{threshold} \
        --mod-threshold 19229:~{threshold} \
        --log-filepath "${filename}.threshold_~{threshold}.tRNA.log" \
        --bedrmod
    >>>

    output {
        File bed = basename(bam, ".bam") + ".threshold_~{threshold}.tRNA.bed"
        File log = basename(bam, ".bam") + ".threshold_~{threshold}.tRNA.log"
    }

    runtime {
        cpu: cpus
        memory: "128GB"
        maxRunTime: 86400 #24 hours (24 * 3600 seconds)
        runtime_minutes: 1 #24 hours (24 * 60 minutes)
        docker: "ontresearch/modkit:sha489d708a48c66368e5d1e118538e5dca68203a64"
        failOnStderr: false
    }
}
