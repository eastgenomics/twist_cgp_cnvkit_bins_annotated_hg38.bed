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
  --refflat       pinned snapshot deposited on DNAnexus (see
                    DEFAULT_REFFLAT_PROJECT_FILE), content-verified against
                    EXPECTED_REFFLAT_SHA256, unless --refflat-path or
                    --refflat-live-ucsc is given

Output:
  twist_cgp_cnvkit_bins_annotated_hg38.bed

Requires: dx (DNAnexus CLI, authenticated), docker (cgp-cnvkit:1.0.0 image
must be docker load'ed from the production image tar -- there is no
registry to pull it from, see the CNVKIT_IMAGE comment below), python3
stdlib only otherwise.
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
# UCSC's live hgdownload URL is not versioned or content-addressed -- it can
# be updated (or removed) at any time with no way to retrieve the exact
# historical file used here. The verified snapshot (see EXPECTED_REFFLAT_SHA256
# below) has been deposited on DNAnexus as an immutable file and is the
# default source; --refflat-live-ucsc is the explicit escape hatch back to
# the live URL (e.g. to deliberately pick up a real UCSC update).
DEFAULT_REFFLAT_PROJECT_FILE = "project-Fkb6Gkj433GVVvj73J7x8KbV:file-JBfX5bQ433Gq7QG1kFPPgpFx"
REFFLAT_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz"

# Content sha256 of the decompressed refFlat.txt, pinned so a future UCSC
# update to the mutable URL above can't silently change the annotation input
# without being noticed -- the build fails loudly instead (see
# --skip-refflat-checksum-assert to intentionally adopt a new refFlat).
# Verified identical across every reproducibility run this file has had
# (original build, two independent from-scratch reruns, and a fresh
# same-day re-download), 2026-09-17.
EXPECTED_REFFLAT_SHA256 = "81815dc035cf03518d144c00f143a36749997dcb5a980c93bf8e346d11e8c800"

# Referenced by immutable content digest (the local image ID from `docker
# load`), not by the mutable tag `cgp-cnvkit:1.0.0` -- a tag can be
# repointed at different image content later without any change here, which
# is exactly the kind of silent drift that produced a different raw
# --annotate output between two runs earlier in this file's history (see
# derive_annotation_normalization_table.py's docstring). This is a local
# image ID (docker inspect --format='{{.Id}}'), not a registry pull digest
# (repo@sha256:...) -- the image is loaded from a tar, not pulled from a
# registry, so there is no registry reference to pin against. `docker run`
# accepts an image ID directly in place of a tag. Confirmed working
# 2026-09-17: `docker run --rm sha256:... cnvkit.py version` -> 0.9.13.
CNVKIT_IMAGE = "sha256:7999995d0d5270ae44e8b130e9cf259e3c60f77831a8bdb46cf2376db24b8aab"

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


def verify_refflat_checksum(path, skip=False):
    import hashlib

    actual = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if skip:
        print(f"refFlat checksum check skipped (--skip-refflat-checksum-assert): {actual}", file=sys.stderr)
        return
    if actual != EXPECTED_REFFLAT_SHA256:
        sys.exit(
            f"refFlat checksum mismatch: expected {EXPECTED_REFFLAT_SHA256}, got {actual}.\n"
            f"UCSC's refFlat.txt (or the file passed via --refflat-path) has changed since "
            f"EXPECTED_REFFLAT_SHA256 was pinned. This would silently change the annotation "
            f"input, so the build stops here instead. If this is an intentional refFlat "
            f"update, review the resulting annotation diff, then update "
            f"EXPECTED_REFFLAT_SHA256 (and re-derive/re-audit the normalization table and "
            f"EXPECTED_BINS_MD5 as needed). To proceed once, without updating the pin, pass "
            f"--skip-refflat-checksum-assert."
        )
    print(f"refFlat checksum verified: {actual}", file=sys.stderr)


def download_refflat(dest, cached_path=None, skip_checksum_assert=False, dnanexus_source=None, use_live_ucsc=False):
    dest = Path(dest)
    if cached_path:
        dest.write_bytes(Path(cached_path).read_bytes())
    elif use_live_ucsc:
        print(
            f"Downloading refFlat.txt.gz from LIVE UCSC (not the pinned DNAnexus "
            f"snapshot -- only the checksum below protects you here): {REFFLAT_URL}",
            file=sys.stderr,
        )
        gz_path = str(dest) + ".gz"
        urllib.request.urlretrieve(REFFLAT_URL, gz_path)
        sh(f"gunzip -kf {gz_path}")
    else:
        source = dnanexus_source or DEFAULT_REFFLAT_PROJECT_FILE
        print(f"Downloading refFlat.txt.gz from the pinned DNAnexus snapshot: {source}", file=sys.stderr)
        gz_path = str(dest) + ".gz"
        dx_download(source, gz_path)
        sh(f"gunzip -kf {gz_path}")
    verify_refflat_checksum(dest, skip=skip_checksum_assert)


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
    comma/semicolon-joined multi-gene labels for overlapping loci.

    Each per-chromosome interval is stored as (start, end, genes,
    prefix_max_end) -- prefix_max_end is the maximum end seen anywhere
    from the start of this chromosome's sorted list up to and including
    this interval. overlapping_genes() uses it to know exactly when a
    backward scan can safely stop; see its docstring."""
    regions = defaultdict(list)
    with open(path) as f:
        for line in f:
            chrom, start, end, name = line.rstrip("\n").split("\t")
            start, end = int(start), int(end)
            base = name.split("_")[0].replace(";", ",")
            genes = set(g for g in base.split(",") if g)
            regions[chrom].append((start, end, genes))
    for c in regions:
        regions[c].sort(key=lambda iv: (iv[0], iv[1]))
        running_max_end = float("-inf")
        with_prefix_max = []
        for s, e, genes in regions[c]:
            running_max_end = max(running_max_end, e)
            with_prefix_max.append((s, e, genes, running_max_end))
        regions[c] = with_prefix_max
    return regions


def overlapping_genes(regions, chrom, start, end):
    """Return the union of gene tokens from every interval on `chrom` that
    overlaps [start, end).

    Intervals are sorted by start, so idx = bisect_left(starts, end) finds
    the point where every interval at index < idx already satisfies
    start < end -- half the overlap condition is guaranteed by position
    alone. The only remaining question for those is whether end > start.

    A fixed look-back window before idx (e.g. checking only the last 5
    entries) can silently miss a real overlap whenever more than that many
    intervals start between the true match and idx -- confirmed against
    this exact panel's targets BED: ALK, MSH2, KIT and PTEN all sit in loci
    dense enough to trigger it with a 5-entry window. Instead, walk
    backward from idx-1 using the precomputed prefix-maximum-end: once the
    maximum end seen from the start of the chromosome's list up to the
    current position is <= start, no interval at or before that position
    can reach past start either, so it's safe to stop.
    """
    import bisect

    ivs = regions.get(chrom, [])
    starts = [iv[0] for iv in ivs]
    idx = bisect.bisect_left(starts, end)
    found = set()

    # Every interval at index < idx already has start < end (guaranteed by
    # bisect_left, since starts is sorted): starts[idx] >= end whenever idx <
    # len(ivs), so nothing from idx onward can ever satisfy start < end. Only
    # the backward scan can find a match.
    k = idx - 1
    while k >= 0:
        s, e, genes, prefix_max_end = ivs[k]
        if e > start:
            found |= genes
        if prefix_max_end <= start:
            break
        k -= 1

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
    ap.add_argument(
        "--refflat-source",
        default=DEFAULT_REFFLAT_PROJECT_FILE,
        help="DNAnexus project:file for the pinned refFlat.txt.gz snapshot",
    )
    ap.add_argument(
        "--refflat-live-ucsc",
        action="store_true",
        help="Download refFlat.txt fresh from live UCSC instead of the pinned DNAnexus "
        "snapshot (still checksum-verified, but UCSC's URL is otherwise unversioned)",
    )
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
    ap.add_argument(
        "--skip-refflat-checksum-assert",
        action="store_true",
        help="Proceed even if refFlat's content doesn't match EXPECTED_REFFLAT_SHA256 "
        "(e.g. while reviewing an intentional UCSC update, before repinning it)",
    )
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
        download_refflat(
            refflat,
            cached_path=args.refflat_path,
            skip_checksum_assert=args.skip_refflat_checksum_assert,
            dnanexus_source=args.refflat_source,
            use_live_ucsc=args.refflat_live_ucsc,
        )
    else:
        for p in (baits_bed, targets_bed, refflat):
            if not p.exists():
                sys.exit(f"--skip-download given but {p} is missing")
        # A reused workdir's refFlat.txt could be stale or hand-edited since it
        # was last downloaded -- --skip-download must not let that reach CNVkit
        # unverified, or the whole point of EXPECTED_REFFLAT_SHA256 is defeated.
        verify_refflat_checksum(refflat, skip=args.skip_refflat_checksum_assert)

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
