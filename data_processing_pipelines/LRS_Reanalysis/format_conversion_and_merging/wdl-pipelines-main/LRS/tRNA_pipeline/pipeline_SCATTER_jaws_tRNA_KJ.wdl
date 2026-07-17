version 1.0

workflow ont_tRNA {
    input {
        Array[File] pod5_files
        String sample_id
        String use_gpu
        Int cpus
        Int gpus
        Boolean barcoded
        String barcoding_model
        String basecalling_type
        String modeldir
    }

    call Seqtagger {
        input:
            pod5_files = pod5_files,
            sample_id = sample_id,
            use_gpu = use_gpu,
            cpus = cpus,
            gpus = gpus,
            barcoded = barcoded,
            barcoding_model = barcoding_model
    }

    scatter (pod5_file in pod5_files) {
        call DoradoBasecall {
            input:
                barcode_table = Seqtagger.barcode_table,
                pod5_file = pod5_file,
                sample_id = sample_id,
                use_gpu = use_gpu,
                cpus = cpus,
                gpus = gpus,
                basecalling_type = basecalling_type,
                modeldir = modeldir
        }
    }

    call MergeBams {
        input:
            bams = DoradoBasecall.bam,
            sample_id = sample_id,
            cpus = cpus
    }

    output {
        File barcode_table = Seqtagger.barcode_table
        Array[File] barcode_pdfs = Seqtagger.pdfs
        Array[File] basecall_bam = DoradoBasecall.bam
        Array[File] barcode_bams = flatten(DoradoBasecall.barcode_bams)
        File merged_bam = MergeBams.merged_bam
        File merged_bai = MergeBams.merged_bai
    }
}


task Seqtagger {
    input {
        Array[File] pod5_files
        String sample_id
        String use_gpu
        Int cpus
        Int gpus
        Boolean barcoded
        String barcoding_model
    }

    command <<<
    set -euo pipefail

    if [ ~{barcoded} = true ]; then
    for pod5 in ~{sep=' ' pod5_files}; do
        filename=$(basename "$pod5" .pod5)
        mkdir -p ./$filename
        mRNA -k /opt/app/models/~{barcoding_model} -r -i $pod5 -o ./$filename/
        zcat ./$filename/*.demux.tsv.gz | awk 'BEGIN { FS = OFS = "\t" } $5 >= 50 { print }' > ./$filename/$filename.tsv
        mv ./$filename/$filename.tsv ./
        mv ./$filename/*.pdf ./
        rm -rf ./$filename
        head -n 1 ./$filename.tsv > header.tsv
    done
    tail -n +2 -q ./*.tsv | sort -k1,1 -k2,2n | cat header.tsv - > "demux_output.tsv"
    else
        echo "Skip" > "demux_output.tsv"
        echo "Skip" > "no_barcoding.pdf"
    fi
    >>>

    output {
        File barcode_table = "demux_output.tsv"
        Array[File] pdfs = glob("*.pdf")
    }

    runtime {
        cpu: cpus
        gpu: true
        gpuCount: gpus
        memory: "1GB"
        maxRunTime: 86400
        runtime_minutes: 1
        docker: "lpryszcz/seqtagger:latest"
    }
}

task DoradoBasecall {
    input {
        File barcode_table
        File pod5_file
        String sample_id
        String use_gpu
        Int cpus
        Int gpus
        String basecalling_type
        String modeldir
    }

    command <<<
    set -euo pipefail
    export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:250
    filename=$(basename "~{pod5_file}" .pod5)
    mkdir -p bam_temp
    mkdir -p final_bams

    if [ ~{basecalling_type} == "sup" ]; then
        dorado basecaller \
            --device ~{use_gpu} \
            "sup,inosine_m6A_2OmeA,m5C_2OmeC,pseU_2OmeU,2OmeG" \
            --models-directory "~{modeldir}" \
            --estimate-poly-a \
            --emit-moves "~{pod5_file}" > "bam_temp/${filename}.bam"
    else
        dorado basecaller \
            --device ~{use_gpu} \
            "fast" \
            --models-directory "~{modeldir}" \
            --estimate-poly-a \
            --emit-moves "~{pod5_file}" > "bam_temp/${filename}.bam"
    fi

    mv "bam_temp/${filename}.bam" "final_bams/${filename}.bam"
    samtools index -b -@ ~{cpus} "final_bams/${filename}.bam"

    if ! grep -q "^Skip$" ~{barcode_table}; then
        mkdir -p "barcode_tables"
        awk -F'\t' -v outdir="barcode_tables" 'NR > 1 {num=sprintf("%02d",$3); print > (outdir "/barcode" num ".tsv")}' ~{barcode_table}
        for tsv in barcode_tables/*.tsv; do
            barcode=$(basename "$tsv" .tsv)
            cut -f 1 "$tsv" > "barcode_tables/${barcode}.read_ids.txt"
            samtools view -@ ~{cpus} -bh -N "barcode_tables/${barcode}.read_ids.txt" "final_bams/${filename}.bam" | samtools sort -@ ~{cpus} -o "final_bams/${barcode}.bam"
            samtools index -b -@ ~{cpus} "final_bams/${barcode}.bam"
        done
        rm -rf "barcode_tables"
    fi

    rm -rf "bam_temp"
    >>>

    output {
        File bam = "final_bams/" + basename(pod5_file, ".pod5") + ".bam"
        Array[File] barcode_bams = glob("final_bams/barcode*.bam")
        File bai = "final_bams/" + basename(pod5_file, ".pod5") + ".bam.bai"
    }

    runtime {
        cpu: cpus
        gpu: true
        gpuCount: gpus
        memory: "1GB"
        maxRunTime: 172800
        runtime_minutes: 1
        docker: "ontresearch/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"
    }
}

task MergeBams {
    input {
        Array[File] bams
        String sample_id
        Int cpus
    }

    command <<<
    set -euo pipefail

    ls ~{sep=' ' bams} > bam_list.txt
    samtools merge -@ ~{cpus} -o "~{sample_id}.merged.bam" -b bam_list.txt
    samtools index -b -@ ~{cpus} "~{sample_id}.merged.bam"
    >>>

    output {
        File merged_bam = "~{sample_id}.merged.bam"
        File merged_bai = "~{sample_id}.merged.bam.bai"
    }

    runtime {
        cpu: cpus
        memory: "128GB"
        maxRunTime: 86400
        runtime_minutes: 1
        docker: "nanozoo/minimap2:2.28--9e3bd01"
    }
}
