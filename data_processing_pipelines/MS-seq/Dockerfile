# Image for running the HRP MS-seq pipelines (NASE_rRNA, NASE_tRNA,
# MS_consensus_creation) on individual input files, e.g. as HPC array jobs.
#
# Build (from this directory):
#   docker build -t hrp-ms-seq:latest .
#
# Base image already ships the OpenMS command line tools required by the
# NASE_rRNA/NASE_tRNA wrappers (DecoyDatabase, NucleicAcidSearchEngine, ...)
# on PATH. See https://github.com/OpenMS/OpenMS/pkgs/container/openms-executables
ARG OPENMS_BASE_IMAGE=ghcr.io/openms/openms-executables:fix_rnome-legacy-bedrmod-export@sha256:0d46a20883df968ed4e1e22e83b7680e708829fe610a10df6dd3b469bf54e9bf
ARG PHASE8_OPENMS_IMAGE=ghcr.io/openms/openms-executables@sha256:f0331d2002d29bc4fc7bc154eb3dc3a3152c8c4b01c4cd37b2474004bda5b7a6
FROM ${PHASE8_OPENMS_IMAGE} AS phase8_decoy
FROM ${OPENMS_BASE_IMAGE}

# Stage 8's RNase-aware DecoyDatabase is required by the archived INIs.
# The newer executable removed the rnase and keep_terminal_nucleotides options.
COPY --from=phase8_decoy /opt/OpenMS/bin/DecoyDatabase /opt/phase8/OpenMS/bin/DecoyDatabase
COPY --from=phase8_decoy /opt/OpenMS/lib/ /opt/phase8/OpenMS/lib/
ENV DECOY_DATABASE_EXECUTABLE=/opt/phase8/OpenMS/bin/DecoyDatabase

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

RUN ln -s /opt/hrp-ms-seq/NASE_rRNA/run_rRNA_analysis.py /usr/local/bin/run_rRNA_analysis \
    && ln -s /opt/hrp-ms-seq/NASE_tRNA/run_tRNA_analysis.py /usr/local/bin/run_tRNA_analysis \
    && ln -s /opt/hrp-ms-seq/MS_consensus_creation/consensus_creation_main.py /usr/local/bin/ms-consensus
    

ENTRYPOINT ["entrypoint.sh"]
CMD ["--help"]
