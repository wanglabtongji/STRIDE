# STRIDE tests

Demo data and an integration test script for STRIDE.

## Data

| File | Dataset | Used by |
| ---- | ------- | ------- |
| `scRNA_E7.5_gene_count.txt`, `scRNA_E7.5_lineage.txt` | Mouse embryo E7.5 scRNA-seq | `embryo` |
| `GEO-seq_E7.5_gene_count.txt`, `GEO-seq_E7.5_location.txt` | Mouse embryo E7.5 GEO-seq (spatial) | `embryo` |
| `Human_heart_scRNA_gene_count.h5`, `Human_heart_scRNA_celltype_curated.txt`, `Human_heart_scRNA_markers.txt` | Developing human heart scRNA-seq | `heart` |
| `Human_heart_ST_gene_count.h5`, `Human_heart_ST_location.txt` | Developing human heart ST | `heart` |

## Running the tests

Install STRIDE (or have it on your `PATH`) and run:

```bash
python test/run_tests.py                 # fast case: mouse embryo (~1-2 min)
python test/run_tests.py --case heart    # human heart: deconvolve + map (~5 min)
python test/run_tests.py --case all      # both
```

Useful flags:

- `--keep` — keep the output files instead of deleting them.
- `--outdir DIR` — write outputs to `DIR` (kept) instead of a temp dir.
- `--stride PATH` — point at a specific `STRIDE` executable.
- `--ntopics N [N ...]` — restrict the topic-number search to speed up a run.

The script exits `0` if all checks pass and `1` otherwise, so it can be wired
into CI.

## What is checked

For each dataset, `STRIDE deconvolve` is run and the
`*_spot_celltype_frac.txt` output is validated:

- the file exists and is non-empty;
- it has spots and cell types, with no `NaN`;
- no spot is left empty (all-zero);
- per-spot fractions sum to 1 (run with `--normalize`);
- all fractions lie in `[0, 1]`.

The empty-`frac`/`NaN` checks specifically guard the regression fixed in
**v1.0.4** (a topic column summing to zero used to propagate `NaN` and yield an
empty result).

For the `heart` case, `STRIDE map` is then run and its outputs are checked to
be present and `NaN`-free, and the mapped cell-type counts are confirmed to
track the deconvolved fractions (consistency fix from v1.0.3/v1.0.4).
