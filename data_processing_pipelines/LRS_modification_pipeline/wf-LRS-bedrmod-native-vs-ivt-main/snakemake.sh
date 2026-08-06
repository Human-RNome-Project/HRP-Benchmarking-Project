#!/bin/bash

module load snakemake/7.32.3 || exit 1

source myconda
snakemake -pr --keep-going --latency-wait 120 --use-conda --conda-frontend conda --use-envmodules -s workflow/Snakefile --profile workflow/profile --rerun-triggers mtime --cluster-cancel scancel