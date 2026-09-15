# Processing of Short-Read Sequencing Data (HRP_B_007_033) Pipeline

# Data Download Instructions
In this step, we will download the raw data files generated directly by the sequencing machines.
1. Open your web browser and navigate to the project's data repository: https://doi.org/10.25585/DOE-HRP/3377574.
2. Follow the instructions on the portal under “Raw data” to download the raw Short-Read Sequencing files (typically .fastq format). Save these files into a dedicated directory on your local machine or high-performance computing (HPC) cluster.

# Requirements & Setup

To run the script you will need a Python 3.8+ environment and R>=4.3.3 with the required dependencies, and an GLORI_pipeline.yml file is provided.

```sh
conda env create -f GLORI_pipeline.yml
conda activate GLORI_pipeline
```

You will also need the GLORI-DUO-tools scrips, which could be download from [GLORI-DUO-tools](https://github.com/ZedekiahZhou/GLORI-DUO-tools)


# Input

- Single-end GLORI *.fastq.gz files with reads which are sense to transcripts
- For paired-end *fastq files, use only the R1 or R2 which is sense to transcripts

# Workflow

1. Prepare analysis pipeline and AG-converted reference for GLORI

First, download the raw reference genome and analysis pipeline

```sh
DUOdir='your_GLORI-DUO-tools_directory'
git clone https://github.com/ZedekiahZhou/GLORI-DUO-tools.git

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

2. Processing raw fastq files, mapping and call m6A using one command

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
  -c 15 -C 5 -s 0.8 -r 0.1 -adp 0.05

# For HRP_B_033, R2 is sense to transcripts; and the 10nt UMI in 5' ends should be clipped
sam="HRP_B_033_2"
rawfq="${sam}_R2.fastq.gz"
outdir="./"
python ${DUOdir}/DUO.py --raw_fq ${rawfq} -o ${outdir} \
  --mode m6A --module preprocessing,mapping,call_m6A --tag_seq none --umi5 10 --nextseq \
  --prx ${sam} \
  -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} \
  -a ${anno} -ba ${baseanno} --gtf ${gtf} --gtf2 ${gtf2} \
  -c 15 -C 5 -s 0.8 -r 0.1 -adp 0.05
```

The output from GLORI-DUO-tools is a txt file `${sam}.totalm6A.FDR.csv.gz` recording detected m6A sites in each row. 

3. Convert output of GLORI-DUO-tools to bedRmod

Use R functions in `convert_to_bedRmod.R` to convert output of GLORI-DUO-tools to bedRmod, and also merge different replicates and methods into the final bedRmod. 


# Output

A merged bedRmod file including m6A modifications from different methods, e.g., GLORI and CAM-seq
