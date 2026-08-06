version 1.0

workflow ont_tRNA {
    input {
        File unaligned_bam
        String sample_id
        File ref_transcriptome
        String basecalling_type
        Int cpus
    }

    call MinimapTranscriptome {
        input:
            bam = unaligned_bam,
            sample_id = sample_id,
            ref_transcriptome = ref_transcriptome,
            cpus = cpus
    }

    if (basecalling_type == "sup" && MinimapTranscriptome.aligned_count > 0) {
        call ModkitPileupTranscriptome {
            input:
                bam = MinimapTranscriptome.aligned_bam,
                bai = MinimapTranscriptome.aligned_bai,
                sample_id = sample_id,
                ref_transcriptome = ref_transcriptome,
                cpus = cpus
        }
    }

    output {
        File transcriptome_bam       = MinimapTranscriptome.aligned_bam
        File transcriptome_bai       = MinimapTranscriptome.aligned_bai
        File? transcriptome_bed      = ModkitPileupTranscriptome.bed
        File? transcriptome_log      = ModkitPileupTranscriptome.log
    }
}


task MinimapTranscriptome {
    input {
        File bam
        String sample_id
        File ref_transcriptome
        Int cpus
    }

    command <<<
    set -euo pipefail

    filename=$(basename "~{bam}" .bam)
    samtools view --threads ~{cpus} -H "~{bam}" | grep "^@RG" > original_rg.txt

    samtools fastq --threads ~{cpus} -T "*" "~{bam}" | \
        minimap2 -y -ax splice -t "~{cpus}" -k 7 -w 3 -A 2 -B 1 -O1,32 -E2,0 -n 1 -m 13 -s 30 --secondary=no --MD "~{ref_transcriptome}" - | \
        samtools sort --threads ~{cpus} | \
        samtools view --threads ~{cpus} -bh -F 260 -o temp.bam

    samtools view --threads ~{cpus} -H temp.bam > new_header.sam
    cat original_rg.txt >> new_header.sam
    samtools reheader new_header.sam temp.bam > "${filename}.transcriptome.aligned.sorted.bam"
    rm temp.bam new_header.sam original_rg.txt
    samtools index -b -@ ~{cpus} "${filename}.transcriptome.aligned.sorted.bam"

    if [ -s "${filename}.transcriptome.aligned.sorted.bam" ]; then
        samtools view --threads ~{cpus} -c "${filename}.transcriptome.aligned.sorted.bam" > aligned_count.txt
    else
        echo 0 > aligned_count.txt
        samtools view -hb -o "${filename}.transcriptome.aligned.sorted.bam" /dev/null
        samtools index -b -@ ~{cpus} "${filename}.transcriptome.aligned.sorted.bam"
    fi
    >>>

    output {
        File aligned_bam   = basename(bam, ".bam") + ".transcriptome.aligned.sorted.bam"
        File aligned_bai   = basename(bam, ".bam") + ".transcriptome.aligned.sorted.bam.bai"
        Int  aligned_count = read_int("aligned_count.txt")
    }

    runtime {
        cpu: cpus
        memory: "128GB"
        maxRunTime: 86400
        runtime_minutes: 1440
        docker: "nanozoo/minimap2:2.28--9e3bd01"
    }
}


task ModkitPileupTranscriptome {
    input {
        File bam
        File bai
        String sample_id
        File ref_transcriptome
        Int cpus
    }

    command <<<
    set -euo pipefail
    filename=$(basename "~{bam}" .bam)
    modkit pileup \
        ~{bam} \
        "${filename}.transcriptome.bed" \
        --ref ~{ref_transcriptome} \
        --threads ~{cpus} \
        --filter-threshold 0.8 \
        --modified-bases m5C 2OmeC inosine m6A 2OmeA pseU 2OmeU 2OmeG \
        --mod-threshold m:0.98 \
        --mod-threshold 19228:0.98 \
        --mod-threshold 17596:0.98 \
        --mod-threshold a:0.98 \
        --mod-threshold 69426:0.98 \
        --mod-threshold 17802:0.98 \
        --mod-threshold 19227:0.98 \
        --mod-threshold 19229:0.98 \
        --log-filepath "${filename}.transcriptome.log" \
        --preload-references \
        --bedrmod
    >>>

    output {
        File bed = basename(bam, ".bam") + ".transcriptome.bed"
        File log = basename(bam, ".bam") + ".transcriptome.log"
    }

    runtime {
        cpu: cpus
        memory: "128GB"
        maxRunTime: 86400
        runtime_minutes: 1440
        docker: "ontresearch/modkit:sha489d708a48c66368e5d1e118538e5dca68203a64"
        failOnStderr: false
    }
}
