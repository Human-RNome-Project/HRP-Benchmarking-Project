# Base image with R 4.2.2 and tidyverse
FROM rocker/tidyverse:4.2.2

LABEL maintainer="Human RNome Project Consortium (HRPC)"
LABEL description="Unified Docker container for MSR-seq tRNA data processing pipeline"
LABEL version="1.0"

ENV DEBIAN_FRONTEND=noninteractive

# 1. Install system dependencies, bioinformatics tools, openjdk for igvtools, and python3
RUN apt-get update && apt-get install -y --no-install-recommends \
    bowtie2 \
    samtools \
    openjdk-17-jre-headless \
    python3 \
    python3-pip \
    python3-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    wget \
    unzip \
    zip \
    && rm -rf /var/lib/apt/lists/*

# 2. Install Python packages (pysam, numpy, pandas)
RUN pip3 install --no-cache-dir \
    pysam \
    numpy \
    pandas

# 3. Install required R packages (openxlsx, BiocManager, Biostrings)
RUN R -e 'install.packages("openxlsx", repos="https://cloud.r-project.org")' && \
    R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install("Biostrings", ask=FALSE)'

# 4. Install IGVTools 2.8.0
RUN wget -q https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_2.8.0.zip -O /tmp/IGV_2.8.0.zip && \
    unzip -q /tmp/IGV_2.8.0.zip -d /opt/ && \
    rm /tmp/IGV_2.8.0.zip && \
    chmod +x /opt/IGV_2.8.0/igvtools && \
    ln -s /opt/IGV_2.8.0/igvtools /usr/local/bin/igvtools

# 5. Copy pipeline scripts into container and place on PATH
RUN mkdir -p /opt/msrseq/scripts
COPY scripts/ /opt/msrseq/scripts/

RUN chmod +x /opt/msrseq/scripts/* && \
    ln -s /opt/msrseq/scripts/refseq_cleanup.R /usr/local/bin/refseq_cleanup.R && \
    ln -s /opt/msrseq/scripts/sam_bin_split.py /usr/local/bin/sam_bin_split.py && \
    ln -s /opt/msrseq/scripts/wig_to_tsv_low_mem_2.py /usr/local/bin/wig_to_tsv_low_mem_2.py && \
    ln -s /opt/msrseq/scripts/tRNAmod_cleaner.R /usr/local/bin/tRNAmod_cleaner.R && \
    ln -s /opt/msrseq/scripts/bedrmod_converter.R /usr/local/bin/bedrmod_converter.R

ENV PATH="/opt/msrseq/scripts:${PATH}"

WORKDIR /data

CMD ["/bin/bash"]
