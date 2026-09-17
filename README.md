# twist_cgp_cnvkit_bins_annotated_hg38.bed

Reproducible build script for `twist_cgp_cnvkit_bins_annotated_hg38.bed` — the
CNVkit target-bin BED for the Twist Oncology DNA Comprehensive Genomic Panel
(CGP) used by `eggd_cgp-cnvkit-{coverage,pon,batch}`.

## What this builds

Twist's capture-probe ("baits") BED is unannotated coordinates only. This
script:

1. Downloads the baits BED and the panel's existing exon-annotated
   ("targets") BED from DNAnexus, and `refFlat.txt` (hg38) from UCSC.
2. Runs `cnvkit.py target --annotate --split` (via the `cgp-cnvkit:1.0.0`
   Docker image already used in production) to add gene names and split
   large bait tiles into properly-sized CNVkit bins.
3. Normalizes CNVkit's raw annotation output against a frozen, audited
   lookup table (`annotation_normalization_table.tsv`) — see
   `derive_annotation_normalization_table.py`'s docstring for why this
   exists (CNVkit's `--annotate` step is not perfectly stable across
   builds).
4. Cross-validates against the targets BED, fills gaps it can, checks
   genuine disagreements against a second independent database (Ensembl),
   and resolves ties by cancer-gene relevance where coordinates alone can't
   decide. Every judgement call is documented inline as an explicit
   `DECISIONS` constant with its evidence — see the module docstring and
   comments in `build_twist_cgp_cnvkit_bins.py`.
5. Leaves anything genuinely unresolved as `UNRESOLVED` rather than
   guessing.

Full narrative (why each step exists, what was tried, what the 131-entry
normalization-table audit found) is in the Confluence controlled document:
[twist_cgp_cnvkit_bins_annotated_hg38.bed](https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/4805623944).

## Usage

```bash
python3 build_twist_cgp_cnvkit_bins.py --out twist_cgp_cnvkit_bins_annotated_hg38.bed
```

Requires: `dx` (authenticated DNAnexus CLI), `docker` (`cgp-cnvkit:1.0.0`
image loaded or pullable), Python 3 stdlib only otherwise.

The build asserts its output against a pinned checksum (`EXPECTED_BINS_MD5`
in the script) — confirmed reproducible across independent from-scratch
runs (fresh downloads, no cached state).

## Tests

```bash
pip install pytest
pytest test_build_twist_cgp_cnvkit_bins.py -v
```

Unit tests cover all pure logic (gene-token parsing, overlap lookups,
normalization-table application, every `DECISIONS` branch, and the
table-derivation script's conflict/alignment checks) with synthetic
fixtures — no Docker, DNAnexus, or network access needed. Integration-level
concerns (`dx_download`, `run_cnvkit_target`, `download_refflat`,
`reverify_ensembl`) are exercised by the reproducibility runs instead, not
mocked here.

## Files

| File | Purpose |
|---|---|
| `build_twist_cgp_cnvkit_bins.py` | Main build script |
| `derive_annotation_normalization_table.py` | Documents/reproduces how the normalization table was derived (diffing a fresh CNVkit run against the original, already-reviewed output) |
| `annotation_normalization_table.tsv` | The frozen, audited ambiguous-name → resolved-name table (131 entries, 14 corrected after an exhaustive cancer-gene-relevance audit — see Confluence doc) |
| `test_build_twist_cgp_cnvkit_bins.py` | Unit tests |
