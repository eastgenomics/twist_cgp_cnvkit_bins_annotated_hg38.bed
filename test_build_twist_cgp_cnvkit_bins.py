#!/usr/bin/env python3
"""
Unit tests for build_twist_cgp_cnvkit_bins.py and
derive_annotation_normalization_table.py.

Covers the pure logic only -- everything that doesn't need Docker, DNAnexus
auth, or network access (dx_download, run_cnvkit_target, download_refflat,
reverify_ensembl are integration-level and exercised by the reproducibility
test instead, see README note below). Every test here uses synthetic BED/TSV
fixtures in a temp directory, so the suite is fast, deterministic, and has
no external dependencies beyond stdlib.

Run with:  python3 -m unittest test_build_twist_cgp_cnvkit_bins.py -v
"""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import build_twist_cgp_cnvkit_bins as build_mod
import derive_annotation_normalization_table as derive_mod


def write_bed(path, rows):
    """rows: list of (chrom, start, end, name)"""
    with open(path, "w") as f:
        for chrom, s, e, name in rows:
            f.write(f"{chrom}\t{s}\t{e}\t{name}\n")


class TestLoadBedWithGenes(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.bed = Path(self.tmp.name) / "targets.bed"

    def tearDown(self):
        self.tmp.cleanup()

    def test_plain_exon_style_name(self):
        write_bed(self.bed, [("chr1", 100, 200, "TNFRSF14_CCDS44046.1_coding_exons_1")])
        regions = build_mod.load_bed_with_genes(self.bed)
        self.assertEqual(regions["chr1"][0][2], {"TNFRSF14"})

    def test_comma_joined_multi_gene_name(self):
        write_bed(self.bed, [("chr13", 100, 200, "ERCC5,BIVM-ERCC5,BIVM")])
        regions = build_mod.load_bed_with_genes(self.bed)
        self.assertEqual(regions["chr13"][0][2], {"ERCC5", "BIVM-ERCC5", "BIVM"})

    def test_semicolon_joined_name_normalized_to_comma(self):
        write_bed(self.bed, [("chr6", 100, 200, "ARID1B;ARID1B")])
        regions = build_mod.load_bed_with_genes(self.bed)
        # both tokens are identical here but the point is ';' must not leak
        # into the gene token itself
        self.assertEqual(regions["chr6"][0][2], {"ARID1B"})

    def test_gene_symbol_with_internal_hyphen_survives_underscore_split(self):
        # H3-3A must not be mistaken for a multi-part underscore name
        write_bed(self.bed, [("chr1", 100, 200, "H3-3A_CCDS1550.1_coding_exons_1")])
        regions = build_mod.load_bed_with_genes(self.bed)
        self.assertEqual(regions["chr1"][0][2], {"H3-3A"})

    def test_rsid_and_clinvar_style_names_kept_as_single_token(self):
        write_bed(
            self.bed,
            [
                ("chr9", 100, 200, "rs1063192"),
                ("chr9", 300, 400, "CDKN2A_RCV000626437"),
            ],
        )
        regions = build_mod.load_bed_with_genes(self.bed)
        self.assertEqual(regions["chr9"][0][2], {"rs1063192"})
        self.assertEqual(regions["chr9"][1][2], {"CDKN2A"})


class TestOverlappingGenes(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.bed = Path(self.tmp.name) / "targets.bed"
        write_bed(
            self.bed,
            [
                ("chr1", 100, 200, "GENEA_exon1"),
                ("chr1", 150, 250, "GENEB_exon1"),  # overlaps GENEA's interval
                ("chr1", 300, 400, "GENEC_exon1"),  # disjoint
            ],
        )
        self.regions = build_mod.load_bed_with_genes(self.bed)

    def tearDown(self):
        self.tmp.cleanup()

    def test_finds_single_overlap(self):
        genes = build_mod.overlapping_genes(self.regions, "chr1", 300, 350)
        self.assertEqual(genes, {"GENEC"})

    def test_finds_multiple_overlapping_records(self):
        genes = build_mod.overlapping_genes(self.regions, "chr1", 160, 190)
        self.assertEqual(genes, {"GENEA", "GENEB"})

    def test_touching_boundary_is_not_an_overlap(self):
        # a bin starting exactly where GENEC's interval ends should not match
        genes = build_mod.overlapping_genes(self.regions, "chr1", 400, 450)
        self.assertEqual(genes, set())

    def test_no_overlap_on_wrong_chromosome(self):
        genes = build_mod.overlapping_genes(self.regions, "chr2", 100, 200)
        self.assertEqual(genes, set())

    def test_no_overlap_in_a_gap(self):
        genes = build_mod.overlapping_genes(self.regions, "chr1", 260, 290)
        self.assertEqual(genes, set())

    def test_finds_overlap_past_more_than_five_intervening_starts(self):
        # Regression test: a fixed 5-entry look-back before the bisect
        # insertion point can silently miss a real overlap when more than
        # five later-starting-but-still-overlapping intervals sit between
        # the true match and the query. This mirrors the real failure
        # found against the actual targets BED (ALK, MSH2, KIT, PTEN).
        # WIDE starts early and runs long; 8 short, densely-packed
        # intervals start after it but well within its span.
        tmp = tempfile.TemporaryDirectory()
        try:
            bed = Path(tmp.name) / "dense.bed"
            rows = [("chr1", 1000, 5000, "WIDE_exon1")]
            for i in range(8):
                s = 2000 + i * 10
                rows.append(("chr1", s, s + 5, f"SHORT{i}_exon1"))
            write_bed(bed, rows)
            regions = build_mod.load_bed_with_genes(bed)
            # Query sits inside WIDE's span but strictly after all 8 SHORT
            # intervals -- so WIDE is 9 positions before the bisect point,
            # well past a 5-entry look-back.
            genes = build_mod.overlapping_genes(regions, "chr1", 4000, 4010)
            self.assertEqual(genes, {"WIDE"})
        finally:
            tmp.cleanup()

    def test_stops_backward_scan_once_prefix_max_end_proves_no_more_overlaps(self):
        # Not a regression test per se -- just confirms the early-exit path
        # itself still returns every genuine overlap, not just the nearest
        # one, once it does apply.
        tmp = tempfile.TemporaryDirectory()
        try:
            bed = Path(tmp.name) / "early_exit.bed"
            write_bed(
                bed,
                [
                    ("chr1", 100, 900, "FAR_exon1"),  # ends well before query
                    ("chr1", 200, 300, "MID_exon1"),
                    ("chr1", 950, 1050, "NEAR_exon1"),  # genuinely overlaps
                ],
            )
            regions = build_mod.load_bed_with_genes(bed)
            genes = build_mod.overlapping_genes(regions, "chr1", 1000, 1100)
            self.assertEqual(genes, {"NEAR"})
        finally:
            tmp.cleanup()


class TestOverlaps(unittest.TestCase):
    def test_genuine_overlap(self):
        self.assertTrue(build_mod.overlaps(("chr1", 100, 200), "chr1", 150, 250))

    def test_touching_is_not_overlap(self):
        self.assertFalse(build_mod.overlaps(("chr1", 100, 200), "chr1", 200, 300))

    def test_different_chromosome_never_overlaps(self):
        self.assertFalse(build_mod.overlaps(("chr1", 100, 200), "chr2", 100, 200))

    def test_fully_contained_interval_overlaps(self):
        self.assertTrue(build_mod.overlaps(("chr1", 100, 200), "chr1", 120, 130))


class TestVerifyRefflatChecksum(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "refFlat.txt"
        self.path.write_bytes(b"some refflat content\n")
        import hashlib

        self.real_sha256 = hashlib.sha256(self.path.read_bytes()).hexdigest()

    def tearDown(self):
        self.tmp.cleanup()

    def test_passes_when_content_matches_pinned_sha256(self):
        orig = build_mod.EXPECTED_REFFLAT_SHA256
        build_mod.EXPECTED_REFFLAT_SHA256 = self.real_sha256
        try:
            build_mod.verify_refflat_checksum(self.path)  # must not raise
        finally:
            build_mod.EXPECTED_REFFLAT_SHA256 = orig

    def test_exits_when_content_does_not_match_pinned_sha256(self):
        orig = build_mod.EXPECTED_REFFLAT_SHA256
        build_mod.EXPECTED_REFFLAT_SHA256 = "0" * 64
        try:
            with self.assertRaises(SystemExit):
                build_mod.verify_refflat_checksum(self.path)
        finally:
            build_mod.EXPECTED_REFFLAT_SHA256 = orig

    def test_skip_flag_bypasses_a_mismatch_without_raising(self):
        orig = build_mod.EXPECTED_REFFLAT_SHA256
        build_mod.EXPECTED_REFFLAT_SHA256 = "0" * 64
        try:
            build_mod.verify_refflat_checksum(self.path, skip=True)  # must not raise
        finally:
            build_mod.EXPECTED_REFFLAT_SHA256 = orig


class TestDownloadRefflatCachedPath(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cached = Path(self.tmp.name) / "cached_refflat.txt"
        self.cached.write_bytes(b"some refflat content\n")
        self.dest = Path(self.tmp.name) / "refFlat.txt"
        import hashlib

        self.real_sha256 = hashlib.sha256(self.cached.read_bytes()).hexdigest()

    def tearDown(self):
        self.tmp.cleanup()

    def test_copies_cached_file_bytes_to_dest(self):
        orig = build_mod.EXPECTED_REFFLAT_SHA256
        build_mod.EXPECTED_REFFLAT_SHA256 = self.real_sha256
        try:
            build_mod.download_refflat(self.dest, cached_path=self.cached)
            self.assertEqual(self.dest.read_bytes(), self.cached.read_bytes())
        finally:
            build_mod.EXPECTED_REFFLAT_SHA256 = orig

    def test_exits_when_cached_file_fails_checksum(self):
        orig = build_mod.EXPECTED_REFFLAT_SHA256
        build_mod.EXPECTED_REFFLAT_SHA256 = "0" * 64
        try:
            with self.assertRaises(SystemExit):
                build_mod.download_refflat(self.dest, cached_path=self.cached)
        finally:
            build_mod.EXPECTED_REFFLAT_SHA256 = orig

    def test_skip_flag_bypasses_a_cached_file_checksum_mismatch(self):
        orig = build_mod.EXPECTED_REFFLAT_SHA256
        build_mod.EXPECTED_REFFLAT_SHA256 = "0" * 64
        try:
            build_mod.download_refflat(
                self.dest, cached_path=self.cached, skip_checksum_assert=True
            )  # must not raise
            self.assertEqual(self.dest.read_bytes(), self.cached.read_bytes())
        finally:
            build_mod.EXPECTED_REFFLAT_SHA256 = orig


class TestDownloadRefflatSourceSelection(unittest.TestCase):
    """download_refflat() must pick DNAnexus (default) vs. live UCSC vs. cached
    correctly, without ever touching the network/dx in a unit test."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dest = Path(self.tmp.name) / "refFlat.txt"
        self.calls = []
        self.orig_dx_download = build_mod.dx_download
        self.orig_urlretrieve = build_mod.urllib.request.urlretrieve
        self.orig_sh = build_mod.sh
        self.orig_sha256 = build_mod.EXPECTED_REFFLAT_SHA256

        def fake_dx_download(source, out_path):
            self.calls.append(("dx_download", source))
            Path(out_path).write_bytes(b"gz-placeholder")

        def fake_urlretrieve(url, out_path):
            self.calls.append(("urlretrieve", url))
            Path(out_path).write_bytes(b"gz-placeholder")

        def fake_sh(cmd, **kw):
            # stands in for `gunzip -kf ...`: just materialise dest with known content
            self.dest.write_bytes(b"some refflat content\n")

        build_mod.dx_download = fake_dx_download
        build_mod.urllib.request.urlretrieve = fake_urlretrieve
        build_mod.sh = fake_sh
        import hashlib

        build_mod.EXPECTED_REFFLAT_SHA256 = hashlib.sha256(b"some refflat content\n").hexdigest()

    def tearDown(self):
        build_mod.dx_download = self.orig_dx_download
        build_mod.urllib.request.urlretrieve = self.orig_urlretrieve
        build_mod.sh = self.orig_sh
        build_mod.EXPECTED_REFFLAT_SHA256 = self.orig_sha256
        self.tmp.cleanup()

    def test_defaults_to_the_pinned_dnanexus_snapshot(self):
        build_mod.download_refflat(self.dest)
        self.assertEqual(self.calls, [("dx_download", build_mod.DEFAULT_REFFLAT_PROJECT_FILE)])

    def test_honours_an_explicit_dnanexus_source_override(self):
        build_mod.download_refflat(self.dest, dnanexus_source="project-X:file-Y")
        self.assertEqual(self.calls, [("dx_download", "project-X:file-Y")])

    def test_live_ucsc_flag_downloads_from_the_ucsc_url_instead(self):
        build_mod.download_refflat(self.dest, use_live_ucsc=True)
        self.assertEqual(self.calls, [("urlretrieve", build_mod.REFFLAT_URL)])

    def test_cached_path_takes_priority_over_both_download_sources(self):
        cached = Path(self.tmp.name) / "cached_refflat.txt"
        cached.write_bytes(b"some refflat content\n")
        build_mod.download_refflat(self.dest, cached_path=cached, use_live_ucsc=True)
        self.assertEqual(self.calls, [])  # neither dx_download nor urlretrieve called


class TestNormalizationTable(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.table_path = Path(self.tmp.name) / "table.tsv"
        self.table_path.write_text(
            "ambiguous_name\tresolved_name\n"
            "MTOR,MTOR-AS1\tMTOR\n"
            "CDKN2B-AS1,CDKN2B\tCDKN2B\n"
        )

    def tearDown(self):
        self.tmp.cleanup()

    def test_loads_table_correctly(self):
        table = build_mod.load_normalization_table(self.table_path)
        self.assertEqual(table["MTOR,MTOR-AS1"], "MTOR")
        self.assertEqual(len(table), 2)

    def test_applies_known_ambiguous_name(self):
        annotated = Path(self.tmp.name) / "annotated.bed"
        write_bed(annotated, [("chr1", 100, 200, "MTOR,MTOR-AS1")])
        table_size, n_unmapped = build_mod.normalize_annotate_output(annotated, self.table_path)
        self.assertEqual(table_size, 2)
        self.assertEqual(n_unmapped, 0)
        result = annotated.read_text().strip().split("\t")
        self.assertEqual(result[-1], "MTOR")

    def test_leaves_clean_names_untouched(self):
        annotated = Path(self.tmp.name) / "annotated.bed"
        write_bed(annotated, [("chr1", 100, 200, "EGFR")])
        build_mod.normalize_annotate_output(annotated, self.table_path)
        result = annotated.read_text().strip().split("\t")
        self.assertEqual(result[-1], "EGFR")

    def test_flags_unmapped_ambiguity_without_crashing(self):
        annotated = Path(self.tmp.name) / "annotated.bed"
        write_bed(annotated, [("chr5", 100, 200, "NEWGENE,NEWGENE-AS1")])
        table_size, n_unmapped = build_mod.normalize_annotate_output(annotated, self.table_path)
        self.assertEqual(n_unmapped, 1)
        # unmapped names are left as-is, not silently guessed
        result = annotated.read_text().strip().split("\t")
        self.assertEqual(result[-1], "NEWGENE,NEWGENE-AS1")


class TestBuildDecisions(unittest.TestCase):
    """Exercises every branch in build_mod.build() against the real
    DECISIONS constants (PROMOTER_CONFIRMED, PSEUDOGENE_OVERRIDE,
    ENSEMBL_H3C14_LOCUS/ENSEMBL_H3C15_LOCUS, DB_DISAGREEMENT_SPANS) so a
    change to any of those is caught by a test, not just a full rebuild."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.targets_bed = Path(self.tmp.name) / "targets.bed"
        self.baits_annotated = Path(self.tmp.name) / "baits_annotated.bed"
        self.out = Path(self.tmp.name) / "out.bed"

    def tearDown(self):
        self.tmp.cleanup()

    def _run(self, annotated_rows, targets_rows):
        write_bed(self.baits_annotated, annotated_rows)
        write_bed(self.targets_bed, targets_rows)
        stats = build_mod.build(self.baits_annotated, self.targets_bed, self.out)
        out_rows = [line.split("\t") for line in self.out.read_text().strip().split("\n")]
        return stats, out_rows

    def test_promoter_confirmed_region_gets_filled(self):
        chrom, s, e = next(iter(build_mod.PROMOTER_CONFIRMED))
        expected_gene = build_mod.PROMOTER_CONFIRMED[(chrom, s, e)]
        stats, rows = self._run(
            annotated_rows=[(chrom, s, e, "-")],
            targets_rows=[],
        )
        self.assertEqual(stats["promoter_confirmed"], 1)
        self.assertEqual(rows[0][3], expected_gene)

    def test_unannotated_region_filled_from_targets_bed(self):
        stats, rows = self._run(
            annotated_rows=[("chr20", 1000, 2000, "-")],
            targets_rows=[("chr20", 1200, 1300, "SOMEGENE_exon1")],
        )
        self.assertEqual(stats["filled_from_targets_bed"], 1)
        self.assertEqual(rows[0][3], "SOMEGENE")

    def test_unannotated_region_with_no_evidence_stays_unresolved(self):
        stats, rows = self._run(
            annotated_rows=[("chr20", 1000, 2000, "-")],
            targets_rows=[],
        )
        self.assertEqual(stats["still_unresolved"], 1)
        self.assertEqual(rows[0][3], "UNRESOLVED")

    def test_pseudogene_override_applied(self):
        ambiguous_name, real_gene = next(iter(build_mod.PSEUDOGENE_OVERRIDE.items()))
        stats, rows = self._run(
            annotated_rows=[("chr1", 226064346, 226064484, ambiguous_name)],
            targets_rows=[],
        )
        self.assertEqual(stats["pseudogene_override"], 1)
        self.assertEqual(rows[0][3], real_gene)

    def test_h3c14_kept_when_inside_real_h3c14_locus(self):
        chrom, s, e = build_mod.ENSEMBL_H3C14_LOCUS
        stats, rows = self._run(
            annotated_rows=[(chrom, s + 10, s + 50, "H3C14")],
            targets_rows=[],
        )
        self.assertEqual(stats["h3_kept_H3C14"], 1)
        self.assertEqual(rows[0][3], "H3C14")

    def test_h3c14_reassigned_to_h3c15_when_inside_real_h3c15_locus(self):
        chrom, s, e = build_mod.ENSEMBL_H3C15_LOCUS
        stats, rows = self._run(
            annotated_rows=[(chrom, s + 10, s + 50, "H3C14")],
            targets_rows=[],
        )
        self.assertEqual(stats["h3_reassigned_to_H3C15"], 1)
        self.assertEqual(rows[0][3], "H3C15")

    def test_db_disagreement_span_overrides_resolved_name(self):
        span, resolved_gene = build_mod.DB_DISAGREEMENT_SPANS[0]
        chrom, s, e = span
        stats, rows = self._run(
            annotated_rows=[(chrom, s + 5, s + 15, "SOME_OTHER_NAME")],
            targets_rows=[],
        )
        self.assertEqual(stats["db_disagreement_resolved"], 1)
        self.assertEqual(rows[0][3], resolved_gene)

    def test_clean_name_outside_all_decisions_passes_through_unchanged(self):
        stats, rows = self._run(
            annotated_rows=[("chr7", 5000, 6000, "EGFR")],
            targets_rows=[],
        )
        self.assertEqual(dict(stats), {})
        self.assertEqual(rows[0][3], "EGFR")

    def test_output_row_count_matches_input(self):
        stats, rows = self._run(
            annotated_rows=[
                ("chr1", 100, 200, "GENEA"),
                ("chr2", 300, 400, "-"),
                ("chr3", 500, 600, "GENEB"),
            ],
            targets_rows=[],
        )
        self.assertEqual(len(rows), 3)


class TestDeriveAnnotationNormalizationTable(unittest.TestCase):
    """derive_annotation_normalization_table.py builds the frozen table by
    diffing a fresh CNVkit run against the original, already-reviewed one --
    these tests cover its own logic in isolation."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.fresh = Path(self.tmp.name) / "fresh.bed"
        self.original = Path(self.tmp.name) / "original.bed"

    def tearDown(self):
        self.tmp.cleanup()

    def test_derives_unambiguous_mapping(self):
        write_bed(self.fresh, [("chr1", 100, 200, "MTOR,MTOR-AS1")])
        write_bed(self.original, [("chr1", 100, 200, "MTOR")])
        table = derive_mod.derive(self.fresh, self.original)
        self.assertEqual(table, {"MTOR,MTOR-AS1": "MTOR"})

    def test_ignores_rows_that_already_match(self):
        write_bed(self.fresh, [("chr1", 100, 200, "EGFR")])
        write_bed(self.original, [("chr1", 100, 200, "EGFR")])
        table = derive_mod.derive(self.fresh, self.original)
        self.assertEqual(table, {})

    def test_ignores_non_ambiguous_disagreements(self):
        # a plain rename with no comma isn't the kind of ambiguity this
        # table is for -- only comma-joined fresh names get recorded
        write_bed(self.fresh, [("chr1", 100, 200, "GENEA")])
        write_bed(self.original, [("chr1", 100, 200, "GENEB")])
        table = derive_mod.derive(self.fresh, self.original)
        self.assertEqual(table, {})

    def test_raises_on_non_unique_resolution(self):
        # the same ambiguous name resolving to two different original
        # values at different positions is a real conflict, not silently
        # resolvable -- must fail loudly, not pick one arbitrarily
        write_bed(
            self.fresh,
            [
                ("chr1", 100, 200, "GENEA,GENEB"),
                ("chr2", 300, 400, "GENEA,GENEB"),
            ],
        )
        write_bed(
            self.original,
            [
                ("chr1", 100, 200, "GENEA"),
                ("chr2", 300, 400, "GENEB"),
            ],
        )
        with self.assertRaises(SystemExit):
            derive_mod.derive(self.fresh, self.original)

    def test_raises_on_position_misalignment(self):
        # fresh and original must be row-aligned (same --split behaviour);
        # a coordinate mismatch means the two files aren't comparable
        write_bed(self.fresh, [("chr1", 100, 200, "GENEA,GENEB")])
        write_bed(self.original, [("chr1", 999, 1000, "GENEA")])
        with self.assertRaises(AssertionError):
            derive_mod.derive(self.fresh, self.original)


if __name__ == "__main__":
    unittest.main()
