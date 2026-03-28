#!/usr/bin/env python3
"""Test suite for i.hyper.spectroscopy — runs without GRASS GIS."""
import sys
import importlib.util
from unittest.mock import MagicMock

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Load the module with grass.script mocked out
# ---------------------------------------------------------------------------

_gs_mock = MagicMock()
sys.modules.setdefault("grass", MagicMock())
sys.modules.setdefault("grass.script", _gs_mock)

_spec = importlib.util.spec_from_file_location(
    "i_hyper_spectroscopy",
    __file__.replace("test_suite.py", "i.hyper.spectroscopy.py"),
)
_mod = importlib.util.module_from_spec(_spec)
sys.modules["i_hyper_spectroscopy"] = _mod  # required for dataclass __module__ lookup
_spec.loader.exec_module(_mod)

# Pull names into module scope for convenience
interpret_spectrum   = _mod.interpret_spectrum
_band_depth          = _mod._band_depth
_compute_depths_numpy   = _mod._compute_depths_numpy
_build_band_species_index = _mod._build_band_species_index
score_composites_numpy   = _mod.score_composites_numpy
SPECTRAL_FEATURES_DB = _mod.SPECTRAL_FEATURES_DB
COMPOSITE_RULES      = _mod.COMPOSITE_RULES
CLASS_LABELS         = _mod.CLASS_LABELS
AbsorptionFeature    = _mod.AbsorptionFeature


# ===========================================================================
# 1. Database integrity
# ===========================================================================

class TestDatabaseIntegrity:
    def test_feature_db_required_fields(self):
        for entry in SPECTRAL_FEATURES_DB:
            assert "center" in entry, f"Missing 'center': {entry}"
            assert "range" in entry, f"Missing 'range': {entry}"
            assert "species" in entry, f"Missing 'species': {entry}"
            assert "mechanism" in entry, f"Missing 'mechanism': {entry}"
            assert "tags" in entry, f"Missing 'tags': {entry}"

    def test_feature_db_wavelength_sanity(self):
        for entry in SPECTRAL_FEATURES_DB:
            lo, hi = entry["range"]
            assert lo < hi, f"Inverted range in {entry['species']}: {lo} >= {hi}"
            assert 380 <= entry["center"] <= 2600, (
                f"Center out of reasonable VNIR+SWIR range: {entry['species']} @ {entry['center']} nm"
            )
            assert lo <= entry["center"] <= hi, (
                f"Center outside its own range: {entry['species']}"
            )

    def test_composite_rules_required_species_exist_in_db(self):
        db_species = {e["species"] for e in SPECTRAL_FEATURES_DB}
        for rule in COMPOSITE_RULES:
            for sp in rule["required_species"]:
                assert sp in db_species, (
                    f"Rule '{rule['name']}': required species '{sp}' not in SPECTRAL_FEATURES_DB"
                )

    def test_composite_rules_bonus_species_exist_in_db(self):
        db_species = {e["species"] for e in SPECTRAL_FEATURES_DB}
        for rule in COMPOSITE_RULES:
            for sp in rule.get("bonus_species", []):
                assert sp in db_species, (
                    f"Rule '{rule['name']}': bonus species '{sp}' not in SPECTRAL_FEATURES_DB"
                )

    def test_composite_rules_base_score_range(self):
        for rule in COMPOSITE_RULES:
            assert 0.0 <= rule["base_score"] <= 1.0, (
                f"Rule '{rule['name']}': base_score {rule['base_score']} out of [0,1]"
            )

    def test_class_labels_match_rules(self):
        assert CLASS_LABELS[0] == "unknown"
        for i, rule in enumerate(COMPOSITE_RULES, start=1):
            assert CLASS_LABELS[i] == rule["name"], (
                f"CLASS_LABELS[{i}] mismatch: '{CLASS_LABELS[i]}' vs '{rule['name']}'"
            )


# ===========================================================================
# 2. _band_depth()
# ===========================================================================

class TestBandDepth:
    def test_simple_dip(self):
        # Feature at 700nm with reflectance 0.2; neighbor at 600nm with 0.8
        depth = _band_depth(0.2, 700.0, [600.0, 700.0, 800.0], [0.8, 0.2, 0.8])
        # continuum = max(0.8, 0.8) = 0.8; depth = 1 - 0.2/0.8 = 0.75
        assert abs(depth - 0.75) < 1e-9

    def test_zero_reflectance(self):
        depth = _band_depth(0.0, 700.0, [600.0, 700.0, 800.0], [0.8, 0.0, 0.8])
        assert abs(depth - 1.0) < 1e-9

    def test_no_neighbors_within_margin(self):
        # All other bands > 150nm away → continuum defaults to 1.0
        depth = _band_depth(0.5, 700.0, [400.0, 700.0, 1000.0], [0.9, 0.5, 0.9])
        # |400-700|=300 > 150, |1000-700|=300 > 150 → no neighbors
        assert abs(depth - 0.5) < 1e-9  # 1 - 0.5/1.0

    def test_feature_equal_to_continuum_gives_zero(self):
        # Feature reflectance matches the continuum → depth = 0
        depth = _band_depth(0.8, 700.0, [600.0, 700.0, 800.0], [0.8, 0.8, 0.8])
        assert depth == 0.0

    def test_depth_not_negative(self):
        # Feature brighter than continuum → clamped to 0
        depth = _band_depth(0.9, 700.0, [600.0, 700.0, 800.0], [0.5, 0.9, 0.5])
        assert depth == 0.0


# ===========================================================================
# 3. interpret_spectrum() — core engine
# ===========================================================================

# Helper to build a flat spectrum with dips
def _make_spectrum(base_wl, base_r, dips):
    """Return triplets list: flat bands at base_r plus dip overrides."""
    points = {wl: base_r for wl in base_wl}
    points.update(dips)
    return [(wl, 20.0, r) for wl, r in sorted(points.items())]


class TestInterpretSpectrumEdgeCases:
    def test_empty_input(self):
        result = interpret_spectrum([])
        assert result.dominant_class is None
        assert result.features == []
        assert result.hypotheses == []

    def test_out_of_range_wavelengths_skipped(self):
        triplets = [(300, 20, 0.5), (2600, 20, 0.5)]  # both out of 400-2500
        result = interpret_spectrum(triplets)
        assert result.dominant_class is None
        assert any("out of 400" in fl for fl in result.flag_notes)

    def test_fwhm_below_threshold_skipped(self):
        # fwhm < 0.5 should be skipped
        triplets = [(700, 0.4, 0.2), (600, 20, 0.8), (800, 20, 0.8)]
        result = interpret_spectrum(triplets)
        # 700nm band is skipped, so no dip → no features
        assert result.dominant_class is None

    def test_all_flat_spectrum_no_features(self):
        triplets = [(wl, 20, 0.5) for wl in range(500, 2400, 100)]
        result = interpret_spectrum(triplets)
        assert result.dominant_class is None
        assert any("No features" in fl for fl in result.flag_notes)

    def test_no_hypothesis_reaches_60_pct(self):
        # Single weak feature that doesn't trigger any composite rule
        triplets = [(700, 20, 0.40), (600, 20, 0.80), (800, 20, 0.80)]
        result = interpret_spectrum(triplets)
        if result.hypotheses and result.hypotheses[0].confidence < 0.60:
            assert result.dominant_class is None
            assert any("60%" in fl for fl in result.flag_notes)

    def test_reflectance_clamped_to_01(self):
        # Reflectance values outside [0,1] should not crash
        triplets = [(600, 20, 1.5), (700, 20, -0.1), (800, 20, 0.8)]
        result = interpret_spectrum(triplets)
        assert result is not None  # no exception

    def test_asymmetry_classification(self):
        # sharp: fwhm < 10
        triplets = [(700, 5, 0.2), (600, 20, 0.8), (800, 20, 0.8)]
        r = interpret_spectrum(triplets)
        if r.features:
            assert r.features[0].asymmetry == "sharp"
        # broad: fwhm > 50
        triplets2 = [(700, 80, 0.2), (600, 20, 0.8), (800, 20, 0.8)]
        r2 = interpret_spectrum(triplets2)
        if r2.features:
            assert r2.features[0].asymmetry == "broad"


class TestInterpretSpectrumMaterials:
    """Physics-based tests: synthetic spectra should classify correctly."""

    def test_green_vegetation(self):
        # chlorophyll_a+b @ 430nm (range 415-445) and 680nm (range 670-695)
        # red_edge_inflection @ 710nm (range 700-730)
        # NIR plateau at 800nm as continuum reference
        triplets = [
            (430, 20, 0.03),   # chlorophyll Soret dip
            (550, 20, 0.15),   # green peak (continuum anchor)
            (680, 20, 0.03),   # chlorophyll red dip
            (710, 20, 0.08),   # red edge dip
            (800, 20, 0.40),   # NIR plateau (continuum anchor)
        ]
        result = interpret_spectrum(triplets)
        assert result.dominant_class == "green_vegetation", (
            f"Expected green_vegetation, got {result.dominant_class}. "
            f"Hypotheses: {[(h.name, h.confidence) for h in result.hypotheses[:5]]}"
        )
        top = result.hypotheses[0]
        assert top.confidence >= 0.60

    def test_liquid_water(self):
        # liquid_water @ 1940nm (range 1900-1980)
        # liquid_water_OH @ 970nm (range 950-990)
        triplets = [
            (800, 20, 0.85),   # baseline NIR
            (970, 20, 0.50),   # liquid_water_OH dip
            (1100, 20, 0.85),  # continuum anchor
            (1800, 20, 0.80),  # baseline SWIR
            (1940, 20, 0.10),  # liquid_water dip (strong)
            (2100, 20, 0.80),  # continuum anchor
        ]
        result = interpret_spectrum(triplets)
        hyp_names = [h.name for h in result.hypotheses]
        assert "liquid_water" in hyp_names, (
            f"Expected liquid_water in hypotheses, got: {hyp_names[:5]}"
        )
        lw = next(h for h in result.hypotheses if h.name == "liquid_water")
        assert lw.confidence >= 0.60

    def test_kaolinite(self):
        # Al-OH_clay @ 1410nm (range 1395-1430)
        # Al-OH @ 2200nm (range 2170-2230)
        triplets = [
            (1300, 20, 0.70),  # continuum
            (1410, 20, 0.45),  # Al-OH_clay dip
            (1560, 20, 0.70),  # continuum
            (2100, 20, 0.70),  # continuum
            (2200, 20, 0.40),  # Al-OH dip
            (2350, 20, 0.70),  # continuum
        ]
        result = interpret_spectrum(triplets)
        hyp_names = [h.name for h in result.hypotheses]
        assert "kaolinite" in hyp_names, (
            f"Expected kaolinite, got: {hyp_names[:5]}"
        )
        k = next(h for h in result.hypotheses if h.name == "kaolinite")
        assert k.confidence >= 0.60

    def test_stressed_vegetation_no_red_edge(self):
        # chlorophyll_a+b present (430nm) but no red_edge_inflection (700-730 absent)
        # → stressed_vegetation, not green_vegetation
        triplets = [
            (430, 20, 0.03),   # chlorophyll Soret dip
            (550, 20, 0.20),   # green peak (continuum)
            (900, 20, 0.40),   # NIR (no 710nm band → no red edge)
        ]
        result = interpret_spectrum(triplets)
        hyp_names = [h.name for h in result.hypotheses]
        # green_vegetation requires red_edge_inflection which is absent
        assert "green_vegetation" not in hyp_names
        assert "stressed_vegetation" in hyp_names

    def test_atmospheric_flag_triggered(self):
        # O2_A_band_atmospheric @ 690nm (range 685-695)
        triplets = [
            (600, 20, 0.80),   # baseline
            (690, 5,  0.50),   # O2 A-band (sharp atmospheric feature)
            (800, 20, 0.80),   # baseline
        ]
        result = interpret_spectrum(triplets)
        assert any("Atmospheric" in fl for fl in result.flag_notes), (
            f"Expected atmospheric flag. Got: {result.flag_notes}"
        )

    def test_hypotheses_sorted_descending(self):
        triplets = [
            (430, 20, 0.03),
            (550, 20, 0.15),
            (680, 20, 0.03),
            (710, 20, 0.08),
            (800, 20, 0.40),
        ]
        result = interpret_spectrum(triplets)
        confs = [h.confidence for h in result.hypotheses]
        assert confs == sorted(confs, reverse=True)

    def test_dominant_class_requires_60pct(self):
        # Green vegetation spectrum → top hypothesis should be ≥ 0.60
        triplets = [
            (430, 20, 0.03),
            (550, 20, 0.15),
            (680, 20, 0.03),
            (710, 20, 0.08),
            (800, 20, 0.40),
        ]
        result = interpret_spectrum(triplets)
        if result.dominant_class is not None:
            assert result.hypotheses[0].confidence >= 0.60


# ===========================================================================
# 4. _compute_depths_numpy()
# ===========================================================================

class TestComputeDepthsNumpy:
    def _make_cube(self, wavelengths, pixel_data):
        """pixel_data: list of reflectance arrays, one per band."""
        return np.array(pixel_data, dtype=np.float64)

    def test_output_shape(self):
        wl = np.array([600., 650., 700., 750., 800.])
        cube = np.ones((5, 10), dtype=np.float64) * 0.5
        depths = _compute_depths_numpy(cube, wl)
        assert depths.shape == (5, 10)

    def test_flat_spectrum_zero_depths(self):
        wl = np.array([600., 650., 700., 750., 800.])
        cube = np.ones((5, 4), dtype=np.float64) * 0.5
        depths = _compute_depths_numpy(cube, wl)
        np.testing.assert_allclose(depths, 0.0, atol=1e-10)

    def test_single_dip_correct_depth(self):
        # 5 bands evenly spaced within 150nm margin of each other
        wl = np.array([600., 650., 700., 750., 800.])
        cube = np.zeros((5, 2), dtype=np.float64)
        cube[:] = 0.80  # baseline
        cube[2, 0] = 0.20  # dip at 700nm in pixel 0
        # For band 2 (700nm), pixel 0:
        # neighbors in (550,850): all other bands qualify
        # continuum = max(0.80, 0.80, 0.80, 0.80) = 0.80
        # depth = 1 - 0.20/0.80 = 0.75
        depths = _compute_depths_numpy(cube, wl)
        assert abs(depths[2, 0] - 0.75) < 1e-9
        # pixel 1 is flat → depth = 0
        assert depths[2, 1] == 0.0

    def test_no_neighbors_within_margin(self):
        # Bands far apart (>150nm) → no neighbors → local_max = 1.0
        wl = np.array([500., 700., 900.])
        cube = np.array([[0.9], [0.5], [0.9]])  # dip at 700nm
        depths = _compute_depths_numpy(cube, wl)
        # 700nm: |700-500|=200 > 150, |700-900|=200 > 150 → no neighbors
        # local_max = 1.0 → depth = 1 - 0.5 = 0.5
        assert abs(depths[1, 0] - 0.5) < 1e-9

    def test_nan_handling(self):
        wl = np.array([600., 650., 700., 750., 800.])
        cube = np.ones((5, 2), dtype=np.float64) * 0.5
        cube[2, 0] = np.nan  # nodata pixel
        depths = _compute_depths_numpy(cube, wl)
        # NaN pixel propagates through the subtraction
        assert np.isnan(depths[2, 0]) or depths[2, 0] >= 0.0  # no crash

    def test_depths_non_negative(self):
        wl = np.linspace(600, 800, 10)
        rng = np.random.default_rng(42)
        cube = rng.uniform(0.0, 1.0, (10, 100))
        depths = _compute_depths_numpy(cube, wl)
        assert np.all(depths >= 0.0)


# ===========================================================================
# 5. _build_band_species_index()
# ===========================================================================

class TestBuildBandSpeciesIndex:
    def test_known_species_found(self):
        # chlorophyll_a+b has ranges (415, 445) and (670, 695)
        wl = np.array([430., 680., 1000.])
        idx = _build_band_species_index(wl)
        assert "chlorophyll_a+b" in idx
        assert 0 in idx["chlorophyll_a+b"]  # 430nm → index 0
        assert 1 in idx["chlorophyll_a+b"]  # 680nm → index 1
        assert 2 not in idx["chlorophyll_a+b"]  # 1000nm out of range

    def test_red_edge_inflection_found(self):
        wl = np.array([710.])  # range (700, 730)
        idx = _build_band_species_index(wl)
        assert "red_edge_inflection" in idx
        assert 0 in idx["red_edge_inflection"]

    def test_no_match_outside_all_ranges(self):
        wl = np.array([350., 2600.])  # outside all DB ranges
        idx = _build_band_species_index(wl)
        # No species should match
        assert len(idx) == 0

    def test_indices_deduplicated_and_sorted(self):
        # Two DB entries for the same species covering overlapping ranges
        wl = np.array([430., 680., 710.])
        idx = _build_band_species_index(wl)
        for sp, band_list in idx.items():
            assert band_list == sorted(set(band_list)), (
                f"Species '{sp}' has duplicate or unsorted indices"
            )


# ===========================================================================
# 6. score_composites_numpy()
# ===========================================================================

class TestScoreCompositesNumpy:
    def test_output_shape(self):
        wl = np.array([430., 550., 680., 710., 800.])
        cube = np.ones((5, 20), dtype=np.float64) * 0.5
        depths = _compute_depths_numpy(cube, wl)
        scores = score_composites_numpy(depths, wl)
        assert scores.shape == (len(COMPOSITE_RULES), 20)

    def test_green_vegetation_pixel_scores(self):
        # pixel 0: dips at 430, 680, 710 → green_vegetation
        # pixel 1: flat → no detection
        wl = np.array([430., 550., 680., 710., 800.])
        cube = np.zeros((5, 2), dtype=np.float64)
        cube[:, 0] = [0.03, 0.15, 0.03, 0.08, 0.40]  # vegetation pixel
        cube[:, 1] = [0.80, 0.80, 0.80, 0.80, 0.80]  # flat pixel
        depths = _compute_depths_numpy(cube, wl)
        scores = score_composites_numpy(depths, wl)

        gv_idx = next(i for i, r in enumerate(COMPOSITE_RULES)
                      if r["name"] == "green_vegetation")
        assert scores[gv_idx, 0] >= 0.60, (
            f"green_vegetation pixel score {scores[gv_idx, 0]:.3f} < 0.60"
        )
        assert scores[gv_idx, 1] == 0.0, (
            f"Flat pixel should score 0.0, got {scores[gv_idx, 1]:.3f}"
        )

    def test_scores_in_range_01(self):
        wl = np.linspace(400, 2500, 50)
        rng = np.random.default_rng(0)
        cube = rng.uniform(0.0, 1.0, (50, 100)).astype(np.float64)
        depths = _compute_depths_numpy(cube, wl)
        scores = score_composites_numpy(depths, wl)
        assert np.all(scores >= 0.0)
        assert np.all(scores <= 1.0)

    def test_required_absent_exclusion(self):
        # stressed_vegetation requires chlorophyll_a+b but NOT red_edge_inflection.
        # If we also provide red_edge_inflection dip → stressed_vegetation should NOT fire.
        wl = np.array([430., 550., 680., 710., 800.])
        cube = np.zeros((5, 1), dtype=np.float64)
        cube[:, 0] = [0.03, 0.15, 0.03, 0.08, 0.40]  # vegetation (has red_edge)
        depths = _compute_depths_numpy(cube, wl)
        scores = score_composites_numpy(depths, wl)
        sv_idx = next(i for i, r in enumerate(COMPOSITE_RULES)
                      if r["name"] == "stressed_vegetation")
        # red_edge_inflection IS detected, so stressed_vegetation is excluded
        assert scores[sv_idx, 0] == 0.0, (
            f"stressed_vegetation should be 0 (red_edge present), got {scores[sv_idx, 0]:.3f}"
        )


# ===========================================================================
# 7. Bonus: feature asymmetry & chromophore scoring
# ===========================================================================

class TestChromophoreScoring:
    def test_confidence_scale_with_depth(self):
        # Deeper absorption → higher chromophore confidence
        shallow = [(680, 20, 0.75), (600, 20, 0.80), (800, 20, 0.80)]
        deep    = [(680, 20, 0.20), (600, 20, 0.80), (800, 20, 0.80)]
        r_shallow = interpret_spectrum(shallow)
        r_deep    = interpret_spectrum(deep)
        # chlorophyll_a (665, range 650-680) should appear with higher confidence for deep
        def _get_chl(result):
            for ca in result.chromophores:
                if "chlorophyll" in ca.species or "hemoglobin" in ca.species or ca.species:
                    return ca
            return None
        # Just check that the deep spectrum has at least as high max confidence
        max_deep    = max((c.confidence for c in r_deep.chromophores), default=0)
        max_shallow = max((c.confidence for c in r_shallow.chromophores), default=0)
        assert max_deep >= max_shallow

    def test_multi_band_bonus_increases_confidence(self):
        # Nd3+ has 4 entries at 580, 745, 800, 867 nm
        # Providing dips at all four should increase chromophore confidence vs one band
        single = [(580, 5, 0.50), (500, 20, 0.80), (650, 20, 0.80)]
        multi  = [
            (500, 20, 0.80), (580, 5, 0.50), (700, 20, 0.80),
            (745, 5, 0.50), (800, 5, 0.50), (867, 5, 0.50),
            (950, 20, 0.80),
        ]
        r_single = interpret_spectrum(single)
        r_multi  = interpret_spectrum(multi)
        def _get_nd(result):
            return next((c for c in result.chromophores if c.species == "Nd3+"), None)
        nd_s = _get_nd(r_single)
        nd_m = _get_nd(r_multi)
        if nd_s and nd_m:
            assert nd_m.confidence >= nd_s.confidence


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
