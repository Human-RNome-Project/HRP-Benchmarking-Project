# I. Input and output

# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz

DUOdir=$1
refdir=$2

raw_genome="${refdir}/GRCh38.primary_assembly.genome.fa.gz"
raw_gtf="${refdir}/gencode.v49.annotation.gtf"
raw_trans="${refdir}/gencode.v49.transcripts.fa.gz"

chr_genome="${refdir}/GRCh38.chr_only.genome.fa"
AGconv_gtf="${refdir}/gencode.v49.annotation.AG_converted.gtf"
longest_trans="${refdir}/gencode.v49.transcripts.longest.fa"


# GLORI reference files

# genome="${refdir}/GRCh38.chr_only.genome.fa.AG_conversion.fa"
# genome2="${refdir}/GRCh38.chr_only.genome.fa"
# rvsgenome="${refdir}/GRCh38.chr_only.genome.fa.rvsCom.fa"
# TfGenome="${refdir}/gencode.v49.transcripts.longest.fa.AG_conversion.fa"
# anno="${refdir}/gencode.v49.annotation.gtf.tbl"
# baseanno="${refdir}/gencode.v49.annotation.gtf.tbl.noredundance.base"
# gtf="${refdir}/gencode.v49.annotation.AG_converted.gtf"
# gtf2="${refdir}/gencode.v49.annotation.gtf"


# I. genome refrence
zcat ${raw_genome} | seqkit grep -r -p "^chr" > ${chr_genome}
python ${DUOdir}/Pre/build_genome_index.py -f ${chr_genome} -p 20 -pre ${chr_genome}


## II. transcriptome refernce
# gtf of only chromosome, gencode v49
cat ${raw_gtf} | awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {$1=$1"_AG_converted"; print}' > ${AGconv_gtf}
python ${DUOdir}/Pre/gtf2anno.py -i ${raw_gtf} -o ${raw_gtf}.tbl

# rename transcripts, replace "|" with " "
zcat ${raw_trans} | sed 's/|/ /g' > ${raw_trans}.renamed
python ${DUOdir}/Pre/selected_longest_transcrpts_fa.py -anno ${raw_gtf}.tbl \
  -fafile ${raw_trans}.renamed --outname_prx ${longest_trans}
python ${DUOdir}/Pre/build_transcriptome_index.py -f ${longest_trans} -pre ${longest_trans}


# III. base anno ======================
python ${DUOdir}/Pre/gtf2genelist.py -i ${raw_gtf} -f ${raw_trans}.renamed -o ${raw_gtf}.genelist
python ${DUOdir}/Pre/anno_to_base.py -i ${raw_gtf}.tbl -o ${raw_gtf}.tbl.baseanno
python ${DUOdir}/Pre/anno_to_base_remove_redundance_v1.0.py -i ${raw_gtf}.tbl.baseanno \
  -o ${raw_gtf}.tbl.noredundance.base -g ${raw_gtf}.genelist
