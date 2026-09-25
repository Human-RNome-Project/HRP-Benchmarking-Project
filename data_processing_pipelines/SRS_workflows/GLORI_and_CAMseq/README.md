# Processing of Short-Read Sequencing Data (HRP_B_007_033) Pipeline

# Data Download Instructions
We will begin by downloading the SRS raw sequencing data based on the metadata available from the [DOE Data Explorer](https://doi.org/10.25585/DOE-HRP/3377574).

1. Open the **Raw Data** tab, followed by the **Metadata** tab, and download `HRP_Metadata_A_LRS.tsv`. In the metadata table, select samples with `mRNA` in the `Sample type` column and either `HRP_B_007` or `HRP_B_033` in the `submission_id` column.

2. The `HRP_B_007` dataset was generated using GLORI. Samples `HRP_B_007_1`, `HRP_B_007_2`, and `HRP_B_007_3` represent three mRNA replicates. Read 1 (R1) corresponds to the sense strand of the transcripts.

3. The `HRP_B_033` dataset was generated using CAM-seq but was analyzed using the GLORI pipeline. Sample `HRP_B_033_1` is an untreated input library and is not used in the current analysis. Samples `HRP_B_033_2` and `HRP_B_033_3` represent two mRNA replicates. Read 2 (R2) corresponds to the sense strand of the transcripts. The 10-nt UMI at the 5′ end of each read must be removed before read mapping.

# Requirements & Setup

To run the script you will need a Python 3.8+ environment and R>=4.3.3 with the required dependencies, and an GLORI_pipeline.yml file is provided.

```sh
conda env create -f GLORI_pipeline.yml
conda activate GLORI_pipeline
```

You will also need the GLORI-DUO-tools scripts, which could be found as `GLORI-DUO-tools.zip`. You could also download it from [GLORI-DUO-tools](https://github.com/ZedekiahZhou/GLORI-DUO-tools)


# Input

- Single-end GLORI *.fastq.gz files with reads which are sense to transcripts
- For paired-end *fastq files, use only the R1 or R2 which is sense to transcripts

# Workflow

1. Step 1: Prepare analysis pipeline and AG-converted reference for GLORI

First, download the raw reference genome and analysis pipeline

```sh
unzip GLORI-DUO-tools.zip
DUOdir='your_GLORI-DUO-tools_directory'

refdir='your_reference_directory'
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
gunzip gencode.v49.annotation.gtf.gz
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz
mv GRCh38.primary_assembly.genome.fa.gz gencode.v49.annotation.gtf gencode.v49.transcripts.fa.gz ${refdir}
```

Then, run the `prepare_reference.sh` to prepare the AG-converted reference. Please make sure `${DUOdir}` refer to your analysis pipeline directory, and the reference file `GRCh38.primary_assembly.genome.fa.gz gencode.v49.annotation.gtf gencode.v49.transcripts.fa.gz` exist in `${refdir}`. 

```sh
bash prepare_reference.sh ${DUOdir} ${refdir}
```

Finally, you will get reference files below:

```sh
genome="${refdir}/GRCh38.chr_only.genome.fa.AG_conversion.fa"
genome2="${refdir}/GRCh38.chr_only.genome.fa"
rvsgenome="${refdir}/GRCh38.chr_only.genome.fa.rvsCom.fa"
TfGenome="${refdir}/gencode.v49.transcripts.longest.fa.AG_conversion.fa"
anno="${refdir}/gencode.v49.annotation.gtf.tbl"
baseanno="${refdir}/gencode.v49.annotation.gtf.tbl.noredundance.base"
gtf="${refdir}/gencode.v49.annotation.AG_converted.gtf"
gtf2="${refdir}/gencode.v49.annotation.gtf"
```

2. Step 2: Processing raw fastq files, mapping and call m6A using one command

```sh
# For HRP_B_007, R1 is sense to trancripts
sam="HRP_B_007_1"
rawfq="${sam}_R1.fastq.gz"
outdir="./"
python ${DUOdir}/DUO.py --raw_fq ${rawfq} -o ${outdir} \
  --mode m6A --module preprocessing,mapping,call_m6A --tag_seq none \
  --prx ${sam} \
  -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} \
  -a ${anno} -ba ${baseanno} --gtf ${gtf} --gtf2 ${gtf2} \
  -c 1 -C 0 -s 0 -r 0 -adp 1.1

# For HRP_B_033, R2 is sense to transcripts; and the 10nt UMI in 5' ends should be clipped
sam="HRP_B_033_2"
rawfq="${sam}_R2.fastq.gz"
outdir="./"
python ${DUOdir}/DUO.py --raw_fq ${rawfq} -o ${outdir} \
  --mode m6A --module preprocessing,mapping,call_m6A --tag_seq none --umi5 10 --nextseq \
  --prx ${sam} \
  -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} \
  -a ${anno} -ba ${baseanno} --gtf ${gtf} --gtf2 ${gtf2} \
  -c 1 -C 0 -s 0 -r 0 -adp 1.1
```

The output from GLORI-DUO-tools is a txt file `${outdir}/03_Sites/${sam}.totalm6A.FDR.csv.gz` recording candidate m6A sites (all A sites with coverage >= 1) in each row. 

3. Step 3: Filter sites and convert output of GLORI-DUO-tools to bedRmod

```sh
# seqMethod: GLORI or CAMseq
Rscript convert_to_bedRmod.R -i ${outdir}/03_Sites/${sam}.totalm6A.FDR.csv.gz \
  -o ${outdir} --seqMethod GLORI -m m6A
```


# Output

A tab-separated bedRMod (`bedRModv2`) with `#key=value` metadata header for a single replicate `${outdir}/${sam}.bed`, with six additonal columns:

- `nRep`: number of replicates that the site is detected, always `1` in the bedRmod for a single replicate
- `repName`: name for each sample (replicate), seperate by ";" for multiple samples
- `repScore`: score for each sample (replicate), seperate by ";" for multiple samples
- `repCov`: coverage for each sample (replicate), seperate by ";" for multiple samples
- `repFreq`: frequency for each sample (replicate), seperate by ";" for multiple samples
- `method`: detection method for each sample (replicate), seperate by ";" for multiple samples