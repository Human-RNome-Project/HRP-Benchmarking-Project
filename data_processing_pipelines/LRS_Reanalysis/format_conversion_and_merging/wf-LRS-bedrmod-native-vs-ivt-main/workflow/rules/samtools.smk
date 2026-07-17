# bgzip_bed coordinate-sorts a bed file (keeping the # header on top) and
# compresses it with BGZF, so the result is tabix-indexable.
rule bgzip_bed:
    input:
        "{filename}.bed"
    output:
        "{filename}.bed.gz"
    threads: 5
    resources:
        mem_mb = 4*1024,
        runtime = 4*60
    conda:
        "../envs/samtools.yml"
    shell:
        """
        tmp="{input}.sorting.tmp"
        {{ grep '^#' {input} || true; }} > "$tmp"
        {{ grep -v '^#' {input} || true; }} \
            | LC_ALL=C sort -k1,1 -k2,2n -T "$(dirname {input})" -S 2G >> "$tmp"
        mv "$tmp" {input}
        bgzip --threads {threads} {input}
        """


# tabix_bed indexes a BGZ-compressed bed file.
rule tabix_bed:
    input:
        "{filename}.bed.gz"
    output:
        "{filename}.bed.gz.tbi"
    threads: 5
    resources:
        mem_mb = 1*1024,
        runtime = 1*60
    conda:
        "../envs/samtools.yml"
    shell:
        """
        tabix --threads {threads} --preset bed {input}
        """
