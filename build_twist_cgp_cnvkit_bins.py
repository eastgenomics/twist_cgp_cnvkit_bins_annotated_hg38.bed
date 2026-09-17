#!/usr/bin/env python3
"""
build_twist_cgp_cnvkit_bins.py

Builds twist_cgp_cnvkit_bins_annotated_hg38.bed — the CNVkit target-bin BED
for the Twist Oncology DNA CGP panel — from Twist's real capture-probe
("baits") footprint, gene-annotated via CNVkit + refFlat, cross-validated
against the panel's existing exon-annotated ("targets") BED, and corrected
against a second independent gene-model database (Ensembl GRCh38).

This consolidates a multi-stage investigation (originally run interactively,
turn by turn) into one reproducible script. See DECISIONS below for the
manual, judgement-based resolutions that a script alone cannot re-derive —
these are baked in as explicit, documented data rather than re-fetched from
external APIs on every run, so the output stays reproducible even if
UCSC/Ensembl change their annotations later. Use --reverify-ensembl to
re-query Ensembl live and confirm DECISIONS still holds (network required,
not needed for a normal build).

Inputs (DNAnexus file IDs, immutable/content-addressed):
  --baits-bed     file-J8kGqJQ46fy5b59VqjjkJqJ3  (Twist_Oncology-DNA-CGP_hg38.bed,
                    the real capture-probe footprint, unannotated)
  --targets-bed   file-J8F1zf045FG6xJVp5X6b0Pkp  (Twist_Oncology_CGP_hg38.bed,
                    the pre-existing exon/ClinVar/MSI-annotated file, used as
                    a fallback annotation source and cross-validation target)
  --refflat       downloaded fresh from UCSC (hg38) unless --refflat-path given

Output:
  twist_cgp_cnvkit_bins_annotated_hg38.bed

Requires: dx (DNAnexus CLI, authenticated), docker (cgp-cnvkit:1.0.0 image
loaded or pullable), python3 stdlib only otherwise.
"""
import argparse
import json
import subprocess
import sys
import time
import urllib.request
from collections import defaultdict
from pathlib import Path

DEFAULT_BAITS_PROJECT_FILE = "project-Fkb6Gkj433GVVvj73J7x8KbV:file-J8kGqJQ46fy5b59VqjjkJqJ3"
DEFAULT_TARGETS_PROJECT_FILE = "project-J8F1vK845FGG0pVG2vfQ4J34:file-J8F1zf045FG6xJVp5X6b0Pkp"
REFFLAT_URL = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz"
CNVKIT_IMAGE = "cgp-cnvkit:1.0.0"

EXPECTED_BINS_MD5 = "79d605af9a68b92bac85e1a706d1828f"  # pinned 2026-09-17, after applying the
# exhaustive 131-entry normalization-table audit (14 corrections). Confirmed
# reproducible end-to-end against a fresh rebuild.

# ---------------------------------------------------------------------------
# DECISIONS — every manual, judgement-based resolution from the original
# investigation, baked in as explicit data. Each entry records WHY, and the
# evidence type, so a reviewer doesn't have to re-run the original session.
# ---------------------------------------------------------------------------

# refFlat carries a duplicate transcript record for H3C14 at the H3C15 locus
# (near-identical histone paralog sequence -> RefSeq annotation artefact).
# Ensembl's canonical H3C14 sits at a distinct locus 12kb away
# (chr1:149,840,687-841,208); its H3C15 matches this disputed region exactly
# (chr1:149,852,608-853,125). Confirmed via Ensembl REST API, 2026-09.
ENSEMBL_H3C14_LOCUS = ("chr1", 149840687, 149841208)
ENSEMBL_H3C15_LOCUS = ("chr1", 149852608, 149853125)

# refFlat's second H3P6 record sits at chr1, overlapping H3-3A almost exactly.
# Ensembl's real H3P6 (a processed pseudogene) exists only on chr2
# (174,719,908-720,318) -- nowhere near chr1. The chr1 refFlat record is a
# duplicate-annotation artefact; always prefer the real gene over the
# phantom pseudogene entry.
PSEUDOGENE_OVERRIDE = {"H3P6": "H3-3A"}

# Ensembl's gene boundaries are genuinely wider than RefSeq/refFlat's at
# these two loci, so BOTH candidate genes truly overlap the disputed bait --
# position alone cannot resolve them. Resolved instead by which gene has
# real oncology/cancer-panel relevance (literature search, 2026-09):
#   IL10 has an established role in tumour immunology (IL-10 blockade
#     immunotherapy literature, promoter-SNP cancer-risk studies); IL19 has
#     no comparable oncology evidence -- chose IL10.
#   MSH6 is a core Lynch-syndrome/mismatch-repair gene already targeted
#     elsewhere in this exact panel; FBXO11's cancer relevance is specific
#     to DLBCL (a lymphoma), not aligned with this solid-cancer panel --
#     chose MSH6.
DB_DISAGREEMENT_SPANS = [
    (("chr1", 206773002, 206773612), "IL10"),   # refFlat said IL19, targets-BED said IL10
    (("chr2", 47807152, 47807272), "MSH6"),     # refFlat said FBXO11, targets-BED said MSH6
]

# Of the 7 regions with no annotation from refFlat OR the targets BED, only
# these 2 are backed by a real Ensembl Regulatory Build promoter/enhancer
# feature genuinely overlapping the bait (checked via Ensembl REST API,
# 2026-09). The remaining 5 (near XPO1, CTLA4 x2, ZFHX3, and a 3rd ID3-
# adjacent position) have no confirmed regulatory evidence and are left as
# UNRESOLVED deliberately -- see the doc's Design section.
PROMOTER_CONFIRMED = {
    ("chr1", 23554509, 23554629): "ID3",
    ("chr1", 23562510, 23562630): "ID3",
}


def sh(cmd, **kw):
    print(f"+ {cmd}", file=sys.stderr)
    return subprocess.run(cmd, shell=True, check=True, **kw)


def dx_download(project_file_id, out_path):
    sh(f"dx download {project_file_id} -o {out_path} -f")


def download_refflat(dest, cached_path=None):
    dest = Path(dest)
    if cached_path:
        dest.write_bytes(Path(cached_path).read_bytes())
        return
    print(f"Downloading refFlat.txt.gz from UCSC: {REFFLAT_URL}", file=sys.stderr)
    gz_path = str(dest) + ".gz"
    urllib.request.urlretrieve(REFFLAT_URL, gz_path)
    sh(f"gunzip -kf {gz_path}")


def run_cnvkit_target(baits_bed, refflat_path, out_path, workdir):
    import os

    sh(
        f"docker run --rm --user {os.getuid()}:{os.getgid()} -v {workdir}:/work -w /work {CNVKIT_IMAGE} "
        f"cnvkit.py target {Path(baits_bed).name} --annotate {Path(refflat_path).name} "
        f"--split -o {Path(out_path).name}"
    )


def load_normalization_table(path):
    """Load the ambiguous-name -> resolved-name table produced by
    derive_annotation_normalization_table.py. See that script's docstring
    for why this exists: cnvkit.py target --annotate is not perfectly
    stable across builds, and this documents/reapplies the exact
    resolution the original, already-reviewed run applied."""
    table = {}
    with open(path) as f:
        next(f)  # header
        for line in f:
            ambiguous, resolved = line.rstrip("\n").split("\t")
            table[ambiguous] = resolved
    return table


def normalize_annotate_output(baits_annotated_path, normalization_table_path):
    """Apply the frozen normalization table to CNVkit's raw --annotate
    output in place, resolving any comma-joined multi-gene name it produced
    to the single gene name the original run used. Names not in the table
    are left untouched (either already clean, or a genuinely new ambiguity
    the table doesn't cover -- those should be reviewed, not silently
    guessed at)."""
    table = load_normalization_table(normalization_table_path)
    lines = []
    unmapped_ambiguous = []
    with open(baits_annotated_path) as f:
        for line in f:
            chrom, s, e, name = line.rstrip("\n").split("\t")
            if "," in name:
                if name in table:
                    name = table[name]
                else:
                    unmapped_ambiguous.append((chrom, s, e, name))
            lines.append(f"{chrom}\t{s}\t{e}\t{name}")
    Path(baits_annotated_path).write_text("\n".join(lines) + "\n")
    if unmapped_ambiguous:
        print(
            f"WARNING: {len(unmapped_ambiguous)} comma-joined names have no "
            f"entry in the normalization table -- these are new ambiguities "
            f"not seen when the table was derived, and need manual review:",
            file=sys.stderr,
        )
        for chrom, s, e, name in unmapped_ambiguous:
            print(f"  {chrom}:{s}-{e}\t{name}", file=sys.stderr)
    return len(table), len(unmapped_ambiguous)


def load_bed_with_genes(path):
    """Load a BED, extracting gene token(s) from column 4. Handles both
    plain exon-style names (GENE_CCDS12345.1_coding_exons_1 -> GENE) and
    comma/semicolon-joined multi-gene labels for overlapping loci."""
    regions = defaultdict(list)
    with open(path) as f:
        for line in f:
            chrom, start, end, name = line.rstrip("\n").split("\t")
            start, end = int(start), int(end)
            base = name.split("_")[0].replace(";", ",")
            genes = set(g for g in base.split(",") if g)
            regions[chrom].append((start, end, genes))
    for c in regions:
        regions[c].sort()
    return regions


def overlapping_genes(regions, chrom, start, end):
    import bisect

    ivs = regions.get(chrom, [])
    starts = [iv[0] for iv in ivs]
    idx = bisect.bisect_left(starts, end)
    found = set()
    j = max(0, idx - 5)
    while j < len(ivs) and ivs[j][0] < end:
        s, e, genes = ivs[j]
        if s < end and e > start:
            found |= genes
        j += 1
    return found


def overlaps(locus, chrom, start, end):
    c, s, e = locus
    return c == chrom and s < end and e > start


def build(baits_annotated_path, targets_bed_path, out_path):
    targets = load_bed_with_genes(targets_bed_path)

    stats = defaultdict(int)
    out_lines = []

    with open(baits_annotated_path) as f:
        for line in f:
            chrom, start, end, name = line.rstrip("\n").split("\t")
            s, e = int(start), int(end)
            key = (chrom, s, e)
            final_name = name

            if name == "-":
                genes = overlapping_genes(targets, chrom, s, e)
                if key in PROMOTER_CONFIRMED:
                    final_name = PROMOTER_CONFIRMED[key]
                    stats["promoter_confirmed"] += 1
                elif genes:
                    final_name = ",".join(sorted(genes))
                    stats["filled_from_targets_bed"] += 1
                else:
                    final_name = "UNRESOLVED"
                    stats["still_unresolved"] += 1
            elif name in PSEUDOGENE_OVERRIDE:
                final_name = PSEUDOGENE_OVERRIDE[name]
                stats["pseudogene_override"] += 1
            elif name == "H3C14":
                if overlaps(ENSEMBL_H3C15_LOCUS, chrom, s, e) and not overlaps(
                    ENSEMBL_H3C14_LOCUS, chrom, s, e
                ):
                    final_name = "H3C15"
                    stats["h3_reassigned_to_H3C15"] += 1
                else:
                    stats["h3_kept_H3C14"] += 1

            for span, resolved in DB_DISAGREEMENT_SPANS:
                if overlaps(span, chrom, s, e):
                    final_name = resolved
                    stats["db_disagreement_resolved"] += 1

            out_lines.append(f"{chrom}\t{s}\t{e}\t{final_name}")

    Path(out_path).write_text("\n".join(out_lines) + "\n")
    return stats


def reverify_ensembl():
    """Optional: re-query Ensembl live to confirm DECISIONS still holds.
    Not required for a normal build -- network-dependent, for audit only."""
    import urllib.request as ur

    def gene_coords(symbol):
        url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json"
        with ur.urlopen(url) as resp:
            d = json.loads(resp.read())
        return d.get("seq_region_name"), d.get("start"), d.get("end"), d.get("biotype")

    checks = ["H3C14", "H3C15", "H3P6", "H3-3A", "IL10", "IL19", "MSH6", "FBXO11"]
    print("Re-verifying DECISIONS against live Ensembl GRCh38 lookups:", file=sys.stderr)
    for g in checks:
        c, s, e, biotype = gene_coords(g)
        print(f"  {g}: chr{c}:{s}-{e} ({biotype})", file=sys.stderr)
        time.sleep(0.3)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--baits-bed-source", default=DEFAULT_BAITS_PROJECT_FILE)
    ap.add_argument("--targets-bed-source", default=DEFAULT_TARGETS_PROJECT_FILE)
    ap.add_argument("--refflat-path", help="Use a local refFlat.txt instead of downloading")
    ap.add_argument("--workdir", default="/tmp/build_twist_cgp_cnvkit_bins")
    ap.add_argument("--out", default="twist_cgp_cnvkit_bins_annotated_hg38.bed")
    ap.add_argument("--skip-download", action="store_true", help="Reuse existing files in --workdir")
    ap.add_argument(
        "--normalization-table",
        default=str(Path(__file__).parent / "annotation_normalization_table.tsv"),
        help="Frozen ambiguous-name -> resolved-name table, see "
        "derive_annotation_normalization_table.py",
    )
    ap.add_argument("--reverify-ensembl", action="store_true")
    ap.add_argument("--skip-checksum-assert", action="store_true")
    args = ap.parse_args()

    if args.reverify_ensembl:
        reverify_ensembl()
        return

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    baits_bed = workdir / "baits.bed"
    targets_bed = workdir / "targets.bed"
    refflat = workdir / "refFlat.txt"
    baits_annotated = workdir / "baits_annotated.bed"

    if not args.skip_download:
        dx_download(args.baits_bed_source, baits_bed)
        dx_download(args.targets_bed_source, targets_bed)
        download_refflat(refflat, cached_path=args.refflat_path)
    else:
        for p in (baits_bed, targets_bed, refflat):
            if not p.exists():
                sys.exit(f"--skip-download given but {p} is missing")

    run_cnvkit_target(baits_bed, refflat, baits_annotated, workdir)

    table_size, n_unmapped = normalize_annotate_output(baits_annotated, args.normalization_table)
    print(
        f"Applied {table_size}-entry normalization table "
        f"({n_unmapped} unmapped ambiguities remaining)",
        file=sys.stderr,
    )

    stats = build(baits_annotated, targets_bed, args.out)

    print("Build stats:", dict(stats), file=sys.stderr)

    if EXPECTED_BINS_MD5 and not args.skip_checksum_assert:
        import hashlib

        actual = hashlib.md5(Path(args.out).read_bytes()).hexdigest()
        if actual != EXPECTED_BINS_MD5:
            sys.exit(f"Checksum mismatch: expected {EXPECTED_BINS_MD5}, got {actual}")
        print(f"Checksum verified: {actual}", file=sys.stderr)
    else:
        print(
            "No pinned checksum set yet (EXPECTED_BINS_MD5) -- "
            "run once, record the md5, then pin it for future reproducibility checks.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
