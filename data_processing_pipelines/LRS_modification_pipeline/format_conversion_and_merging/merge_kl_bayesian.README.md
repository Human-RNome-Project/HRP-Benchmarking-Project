# merge_kl_bayesian.sh

Left-joins the scores from a **Bayesian** and a **KL** bedRMod file onto a **native** bedRMod file. The result is the native file unchanged, with two extra columns at the end: `score_bayesian` and `score_kl`.

The script uses an external sort followed by `join`, so memory use stays at about 2 GB no matter how large the input files are. It handles files with billions of rows.

## Usage

```bash
bash merge_kl_bayesian.sh BAYES_FILE KL_FILE NATIVE_FILE OUTPUT_FILE
```

The order of the arguments matters:

| # | Argument | Description |
|---|----------|-------------|
| 1 | `BAYES_FILE` | bedRMod file with the Bayesian scores (lookup) |
| 2 | `KL_FILE` | bedRMod file with the KL scores (lookup) |
| 3 | `NATIVE_FILE` | native file; every row is kept (left side of the join) |
| 4 | `OUTPUT_FILE` | file to write; must not be one of the inputs |

Example:

```bash
bash merge_kl_bayesian.sh \
  all_mods_pileup_rRNA_native_bayes_modscore.bedrmod \
  all_mods_pileup_rRNA_native_kl_modscore.bedrmod \
  rRNA.native_vs_ivt_fisher.bed \
  rRNA_native_bayesian_KL_joined.bedrmod
```

### Optional environment variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `SORT_MEM` | `2G` | Memory buffer for `sort`. More is faster if the node has the RAM. |
| `SCRATCH` | current directory | Where temporary files are written. |

```bash
SORT_MEM=8G bash merge_kl_bayesian.sh bayes.bedrmod kl.bedrmod native.bed out.bedrmod
```

## Requirements

- bash, plus GNU `awk`, `sort`, `join` and `cut` (standard on Linux clusters)
- Free disk space in `$SCRATCH` of about the combined size of the input files, for the temporary sorted files. They are deleted automatically when the script finishes, including if it fails.

## Input format

All three files are tab-separated bedRMod-style files. Lines starting with `#` are treated as headers.

**Join key.** Rows are matched on five columns:

| Column | Field |
|--------|-------|
| 1 | chrom |
| 2 | chromStart |
| 3 | chromEnd |
| 4 | name (modification) |
| 6 | strand |

Matches must be exact, so all files need the same coordinate convention (0-based starts) and the same chromosome names.

**Modification names.** The Bayesian and KL files may use modkit codes in the `name` column (e.g. `a`, `m`, `17596`). The script translates them to readable names using the `#modification_names=` header line of each file. For example, the header

```
#modification_names=m:m5C:C,19228:Cm:C,17596:I:A,a:m6A:A,...
```

translates `a` → `m6A`, `17596` → `I`, and so on. Codes not listed in the header are left as they are. The native file's `name` column must use the readable names (`m6A`, `Gm`, `Y`, …).

**Scores.** Column 5 (`score`) is taken from the Bayesian and KL files.

## Output

- All header lines of the native file are copied.
- In the column-name header line (`#['chrom', 'chromStart', ...]`), `'score'` is renamed to `'score_native'`, and `'score_bayesian'` and `'score_kl'` are appended at the end of the list.
- Every native row appears once, in its original order, with two extra columns: the Bayesian score, then the KL score.
- A `.` means no matching row was found in that lookup file.

Note that Bayesian and KL scores are on the BED scale (0–1000), while the native score is a fraction (0–1).

## Checking the result

```bash
# row count: must equal the native file
grep -vc '^#' NATIVE_FILE OUTPUT_FILE

# how many rows got a score from each lookup file
awk -F'\t' '!/^#/ { n++; if ($(NF-1)!=".") b++; if ($NF!=".") k++ }
  END { print "rows:", n; print "with Bayesian score:", b+0; print "with KL score:", k+0 }' OUTPUT_FILE
```

- **More rows than the native file** means a key occurs more than once in a lookup file. `join` outputs one row per match, so duplicates multiply rows. Check for duplicate keys with:

  ```bash
  awk -F'\t' '!/^#/ {print $1"|"$2"|"$3"|"$4"|"$6}' BAYES_FILE | LC_ALL=C sort -S 2G | uniq -d | head
  ```

  This shows the raw codes; translate them using the header if needed.

- **Zero or very few matches** usually means a key column differs between the files: a different `name` spelling, chromosome naming (`chr1` vs `1`), or a coordinate offset (1-based vs 0-based). Compare the first rows of each file with:

  ```bash
  grep -v '^#' FILE | head -3 | cut -f1-6 | cat -A
  ```

## How it works

1. Each lookup file is reduced to `key <TAB> score` and sorted by key.
2. The native file is reduced to `key <TAB> line number <TAB> full line` and sorted by key.
3. Two `join -a 1` passes attach the Bayesian and KL scores. `-a 1` keeps unmatched native rows; `-e .` fills missing values with `.`.
4. The result is re-sorted by the saved line number to restore the original order, and the helper columns are removed.

`sort` works externally: it sorts chunks that fit into `SORT_MEM`, writes them to disk and merges them. `join` reads both sorted inputs line by line. Neither ever holds a whole file in memory. This trades RAM for disk space and runtime compared to an in-memory (hash) join.

`LC_ALL=C` is set so that `sort` and `join` use the same byte-wise ordering. Without it, locale-dependent ordering can make `join` silently miss matches.
