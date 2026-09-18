# twist_cgp_cnvkit_bins_annotated_hg38.bed

Reproducible build script for `twist_cgp_cnvkit_bins_annotated_hg38.bed` — the
CNVkit target-bin BED for the Twist Oncology DNA Comprehensive Genomic Panel
(CGP) used by `eggd_cgp-cnvkit-{coverage,pon,batch}`.

## What this builds

Twist's capture-probe ("baits") BED is unannotated coordinates only. This
script:

1. Downloads the baits BED and the panel's existing exon-annotated
   ("targets") BED from DNAnexus, plus `refFlat.txt` (hg38) via one of four
   source modes (see the table below) -- UCSC's live hgdownload file is
   unversioned and could change or disappear without notice, so the exact
   historical snapshot used here is deposited on DNAnexus and used by
   default.
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

### refFlat source modes

Checked in this order — the first that applies wins; every mode still goes
through the same `EXPECTED_REFFLAT_SHA256` content check unless
`--skip-refflat-checksum-assert` is also passed:

| Precedence | Option | Behaviour |
|---|---|---|
| 1 (highest) | `--refflat-path PATH` | Use a local file directly, no download |
| 2 | `--refflat-live-ucsc` | Download fresh from the live UCSC URL (`REFFLAT_URL`) |
| 3 | `--refflat-source PROJECT:FILE` | Download this specific DNAnexus file instead of the default |
| 4 (default) | *(none of the above given)* | Download the pinned DNAnexus snapshot, `DEFAULT_REFFLAT_PROJECT_FILE` |

`--refflat-source` and `--refflat-live-ucsc` are mutually exclusive in
practice: `--refflat-live-ucsc` takes priority whenever both are given,
since it's checked first.

## Usage — recreating the file

```bash
python3 build_twist_cgp_cnvkit_bins.py --out twist_cgp_cnvkit_bins_annotated_hg38.bed
```

Requires: `dx` (authenticated DNAnexus CLI), `docker` (`cgp-cnvkit:1.0.0`
image must be `docker load`ed from the production image tar -- there is no
registry to pull it from), Python 3 stdlib only otherwise.

Run in a fresh/empty `--workdir` (the default,
`/tmp/build_twist_cgp_cnvkit_bins`, is fine as long as nothing from a
previous run is left in it), without `--skip-download`, `--refflat-path`,
or any other cached-input override. This is the exact same command whether
you're recreating the file for the first time or independently verifying
reproducibility -- there is no separate test mode or test script. A real,
fresh run downloads the baits BED, targets BED, and refFlat from scratch,
runs them through `cnvkit.py target --annotate --split`, and the script
asserts the result against the pinned checksum below (`EXPECTED_BINS_MD5`),
exiting non-zero on any mismatch -- so a clean exit *is* the reproducibility
proof. Confirmed this way across multiple independent from-scratch runs
during development.

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
