# Image for running the HRP MS-seq pipelines (NASE_rRNA, NASE_tRNA,
# MS_consensus_creation) on individual input files, e.g. as HPC array jobs.
#
# Build (from this directory):
#   docker build -t hrp-ms-seq:latest .
#
# Base image already ships the OpenMS command line tools required by the
# NASE_rRNA/NASE_tRNA wrappers (DecoyDatabase, NucleicAcidSearchEngine, ...)
# on PATH. See https://github.com/OpenMS/OpenMS/pkgs/container/openms-executables
FROM ghcr.io/openms/openms-executables:RNOME

USER root

# consensus_creation_main.py needs pandas (see MS_consensus_creation/environment.yml).
# The NASE_rRNA/NASE_tRNA wrapper scripts only need the python3 stdlib plus the
# OpenMS executables already provided by the base image.
RUN apt-get update \
    && apt-get install -y --no-install-recommends python3 python3-venv \
    && rm -rf /var/lib/apt/lists/* \
    && python3 -m venv /opt/venv

ENV PATH="/opt/venv/bin:${PATH}"

RUN pip install --no-cache-dir pandas

WORKDIR /opt/hrp-ms-seq
COPY MS_consensus_creation/*.py ./MS_consensus_creation/
COPY NASE_rRNA/ ./NASE_rRNA/
COPY NASE_tRNA/ ./NASE_tRNA/
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh

RUN chmod +x /usr/local/bin/entrypoint.sh \
    NASE_rRNA/run_rRNA_analysis.py \
    NASE_tRNA/run_tRNA_analysis.py \
    MS_consensus_creation/consensus_creation_main.py

ENTRYPOINT ["entrypoint.sh"]
CMD ["--help"]
