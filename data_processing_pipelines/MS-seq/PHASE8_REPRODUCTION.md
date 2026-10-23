# Reproduce the saved phase 8 MS-seq searches

For the tRNA phase 8 reference, all three samples used the same RNase-aware RNA decoy database; the earlier Kaiser source INI had `type=protein` and does not reproduce that reference.

The `NASE_rRNA/phase8_profiles` and `NASE_tRNA/phase8_profiles` directories contain the validated, sample-specific NucleicAcidSearchEngine and DecoyDatabase INIs. The default INIs in `NASE_rRNA` and `NASE_tRNA` are the Dedon profiles. Use the profile for the corresponding sample; the labs did not all use the same search settings.

Build the Dockerfile in this directory. It pins the newer NucleicAcidSearchEngine with the restored observation-level bedRMod writer and the phase 8 DecoyDatabase executable, which supports the RNase-aware decoy settings in these INIs. Both source images are pinned by digest. The INIs request 28 search threads.

For an exact phase 8 reproduction, pass the sample's profile with `--nase-ini` and `--decoy-ini`. The search tolerances are already in its NASE.ini, so omit the command-line tolerance flags. If the flags are supplied, their values take precedence exactly as supplied.

| Sample | Mode | Precursor ppm in phase 8 | Fragment ppm in phase 8 |
| --- | --- | ---: | ---: |
| HRP_C_001_1 | rRNA | 15 | 14 |
| HRP_C_001_3 | tRNA | 10 | 10 |
| HRP_C_002_1 | rRNA | 5 | 5.5 |
| HRP_C_005_4 | rRNA | 10 | 30 |
| HRP_C_006_4 | rRNA | 10 | 30 |
| HRP_C_008_4 | tRNA | 10 | 15 |
| HRP_C_014_1 | rRNA | 20 | 20 |
| HRP_C_015_2 | tRNA | 20 | 20 |

For example, for `HRP_C_001_1` with the project mounted at `/data`:

```sh
docker run --rm -v "$PWD":/data hrp-ms-seq:latest rrna \
  /data/Mzml_data/HRP_C_001_1.mzML /data/hs_rRNA_Us_28S+18S.fa \
  --output-dir /data/Mzml_data/phase8_results \
  --nase-ini /opt/hrp-ms-seq/NASE_rRNA/phase8_profiles/HRP_C_001_1/NASE.ini \
  --decoy-ini /opt/hrp-ms-seq/NASE_rRNA/phase8_profiles/HRP_C_001_1/decoy_database.ini
```

The tolerances in the previously supplied grid-search table differ from phase 8 for six of the eight samples. Running those explicit values is valid for a new grid search, but it cannot be expected to reproduce the saved phase 8 files byte for byte. The CSV mapping bundled for tRNA also matches the phase 8 mapping (`mchm5U` has ID `99998`). Exact reproduction additionally requires the same mzML and target FASTA content.
