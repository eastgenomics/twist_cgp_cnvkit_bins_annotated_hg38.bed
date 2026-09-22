#!/usr/bin/env python3
"""
derive_annotation_normalization_table.py

Documents and reproduces HOW annotation_normalization_table.tsv was derived.

CNVkit's `cnvkit.py target --annotate` step is not perfectly stable across
different installations/builds of the same reported version (0.9.13): a
fresh run against byte-identical inputs (same baits BED, same refFlat.txt,
same "cgp-cnvkit:1.0.0" image tag) produces comma-joined multi-gene names
for 662 of 14,393 regions where the *original* run -- whose output is the
one already cross-validated, Ensembl-corrected, and published as the
approved twist_cgp_cnvkit_bins_annotated_hg38.bed -- produced a single clean
gene name. Root cause not fully pinned down (ruled out: input differences,
PYTHONHASHSEED); most likely an undocumented drift in the exact CNVkit git
commit baked into the Docker image between the two builds.

Rather than treat this as noise or write a hand-guessed heuristic rule (a
"drop antisense/LOC/MIR names" rule does NOT hold universally -- see
CDKN2B-AS1,CDKN2B -> CDKN2B-AS1, where the antisense transcript name is the
one that was originally kept), this script mechanically extracts the exact
mapping the original run already applied, by diffing the two raw annotate
outputs line-for-line. This makes the correction data-derived and
reviewable rather than an inferred rule that could silently mis-resolve an
edge case.

Usage:
    python3 derive_annotation_normalization_table.py \
        --fresh baits_annotated_fresh.bed \
        --original baits_annotated_original.bed \
        --out annotation_normalization_table.tsv

Re-run this only if the CNVkit build changes again and a new drift needs
documenting -- do not re-run casually, since it hard-codes whatever the
"original" file says, which is only valid because the original file is the
already-reviewed baseline.
"""
import argparse
import csv
from collections import defaultdict
from itertools import zip_longest

_MISSING = object()  # sentinel: a row present on one side but not the other


def derive(fresh_path, original_path, stats=None):
    """Derive the ambiguous-name -> resolved-name table.

    `stats`, if given a dict, is populated in-place with `n_occurrences`:
    the total number of individual regions that hit an ambiguous name
    (as opposed to the number of *distinct* ambiguous name strings, which
    is just `len()` of the returned table) -- e.g. 662 regions collapsing
    down to 131 unique ambiguous names.
    """
    mapping = defaultdict(set)  # ambiguous_name -> {resolved single names}
    n_occurrences = 0
    with open(fresh_path) as ff, open(original_path) as fo:
        for i, (fresh_line, orig_line) in enumerate(zip_longest(ff, fo, fillvalue=_MISSING), start=1):
            # zip_longest (not zip) so a line-count mismatch between the two
            # files raises loudly here instead of zip silently truncating to
            # the shorter file and dropping trailing rows unnoticed.
            if fresh_line is _MISSING or orig_line is _MISSING:
                shorter, longer = ("fresh", "original") if fresh_line is _MISSING else ("original", "fresh")
                raise SystemExit(
                    f"Line count mismatch -- fresh and original files must have the "
                    f"same number of rows (same baits BED, same --split behaviour): "
                    f"{shorter} file has exactly {i - 1} lines, {longer} file has at least {i} lines"
                )
            f_chrom, f_s, f_e, f_name = fresh_line.rstrip("\n").split("\t")
            o_chrom, o_s, o_e, o_name = orig_line.rstrip("\n").split("\t")
            assert (f_chrom, f_s, f_e) == (o_chrom, o_s, o_e), (
                f"Position mismatch -- fresh and original files must be "
                f"row-aligned (same baits BED, same --split behaviour): "
                f"{f_chrom}:{f_s}-{f_e} vs {o_chrom}:{o_s}-{o_e}"
            )
            if f_name != o_name and "," in f_name:
                mapping[f_name].add(o_name)
                n_occurrences += 1
    # every ambiguous name must resolve to exactly one original value --
    # if it doesn't, the mapping is ambiguous itself and needs a human look
    conflicts = {k: v for k, v in mapping.items() if len(v) > 1}
    if conflicts:
        raise SystemExit(f"Non-unique resolutions found, needs manual review: {conflicts}")
    if stats is not None:
        stats["n_occurrences"] = n_occurrences
    return {k: next(iter(v)) for k, v in mapping.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fresh", required=True)
    ap.add_argument("--original", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    stats = {}
    table = derive(args.fresh, args.original, stats=stats)

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["ambiguous_name", "resolved_name"])
        for ambiguous, resolved in sorted(table.items()):
            w.writerow([ambiguous, resolved])

    print(
        f"Found {stats['n_occurrences']} ambiguous regions, "
        f"resolving to {len(table)} unique normalization entries -- wrote to {args.out}"
    )


if __name__ == "__main__":
    main()
