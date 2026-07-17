# wf-ont-bedrmod-native-vs-ivt

Snakemake workflow for identifying RNA modifications from Oxford Nanopore
direct-RNA sequencing data by comparing native and in-vitro transcribed (IVT)
samples.

---

## Methods

### Identification of RNA modification sites

To distinguish genuine RNA modifications from sequencing noise, modification
calls derived from native direct-RNA sequencing were compared against calls from
an in-vitro transcribed (IVT) control, which is unmodified by design. For each
genomic position present in the native BEDRMod file, a one-sided Fisher's exact
test was performed to assess whether the modification rate in the native sample
was significantly greater than in the IVT control. Specifically, for a given
site, a 2×2 contingency table was constructed from the counts of modified and
canonical reads in each condition, and the p-value was computed as P(X ≥
n<sub>mod,nat</sub>), where X follows a Hypergeometric(N, K, n) distribution
with N = total reads across both conditions, K = total modified reads, and n =
total native reads. Only sites with a minimum coverage of 30 reads in both the
native and IVT samples were considered testable; sites below this threshold were
retained in the output but assigned NA for all statistical fields. Sites absent
from the IVT sample were similarly assigned NA p-values. Multiple testing
correction was applied to all testable sites jointly using the
Benjamini–Hochberg false discovery rate (FDR) procedure.
