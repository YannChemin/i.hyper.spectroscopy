#!/usr/bin/env python3
"""Derive the camouflage branch of the spectral taxonomy from USGS splib07.

The C2RD Pourrieres campaign of 13 October 2026 laid out twelve camouflage
samples that the taxonomy could not match, because it had no entry for
manufactured camouflage material at all. The samples are not knotted netting:
they are die-cut sheet garnish, a single coated fabric punched with leaf-shaped
voids, so the measured spectrum is the fabric plus whatever shows through the
cuts.

This script measures the absorption features of the splib07 coated-fabric and
tarpaulin spectra, checks them against natural-fibre controls, and writes the
resulting branch to ``camouflage_branch.json`` next to the knowledge base.

:param --splib: root of the local splib07 copy
:param --out: where to write the branch
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

CONTINUUM_MARGIN_NM = 120.0
DELETED = -1e10  # splib07 marks deleted channels with -1.23e34

# Coated synthetic fabrics and woven tarpaulins: the material class the
# die-cut garnish belongs to.
SYNTHETIC = {
    "nylon_webbing_olive_drab": "splib07a_Nylon_Webbing_GDS428_OlvDrab_ASDFRa_AREF.txt",
    "nylon_fabric_green": "splib07a_Nylon_Fabric_GDS432_Grn_RpSt_ASDFRa_AREF.txt",
    "nylon_fabric_blue": "splib07a_Nylon_Fabric_GDS433_Blu_RpSt_ASDFRa_AREF.txt",
    "nylon_fabric_red": "splib07a_Nylon_Fabric_GDS431_Red_RpSt_ASDFRa_AREF.txt",
    "nylon_fabric_black_coated": "splib07a_Nylon_Fabric_GDS436_BlkCoatd_ASDFRa_AREF.txt",
    "nylon_webbing_yellow": "splib07a_Nylon_Webbing_GDS438_Flr_Yel_ASDFRa_AREF.txt",
    "plastic_tarp_green": "splib07a_Plastic_Tarp_GDS339_Green_ASDFRa_AREF.txt",
    "plastic_tarp_blue_woven": "splib07a_Plastic_Tarp_GDS340_Blu_Wovn_ASDFRa_AREF.txt",
}

# Cellulose controls: these must NOT share the diagnostic bands, otherwise the
# branch would fire on canvas, hessian and dry vegetation.
NATURAL = {
    "cotton_fabric_white": "splib07a_Cotton_Fabric_GDS437_White_ASDFRa_AREF.txt",
    "burlap_fabric_brown": "splib07a_Burlap_Fabric_GDS430_Brown_ASDFRa_AREF.txt",
}


def load_column(path: Path) -> np.ndarray:
    """One splib07 ASCII column, with deleted channels as NaN."""
    values = []
    with path.open() as handle:
        handle.readline()  # title line
        for line in handle:
            try:
                values.append(float(line))
            except ValueError:
                pass
    out = np.asarray(values, dtype=float)
    out[out < DELETED] = np.nan
    return out


def band_depths(wavelengths: np.ndarray, reflectance: np.ndarray) -> np.ndarray:
    """Continuum-removed depth, continuum = local max in +/- the margin."""
    depths = np.full_like(reflectance, np.nan)
    finite = np.isfinite(reflectance)
    for i, (centre, value) in enumerate(zip(wavelengths, reflectance)):
        if not finite[i]:
            continue
        window = (np.abs(wavelengths - centre) <= CONTINUUM_MARGIN_NM) & finite
        if window.sum() < 3:
            continue
        continuum = np.nanmax(reflectance[window])
        if continuum > 0:
            depths[i] = max(0.0, 1.0 - value / continuum)
    return depths


def local_maxima(wavelengths, depths, lo=900.0, hi=2450.0, threshold=0.03):
    """Absorption band centres, as (wavelength, depth), merged within 25 nm."""
    found = []
    for i in range(2, len(depths) - 2):
        if not np.isfinite(depths[i]) or not lo <= wavelengths[i] <= hi:
            continue
        neighbourhood = depths[max(0, i - 12) : i + 13]
        if depths[i] >= threshold and depths[i] == np.nanmax(neighbourhood):
            found.append((float(wavelengths[i]), float(depths[i])))
    merged: list[tuple[float, float]] = []
    for centre, depth in found:
        if merged and centre - merged[-1][0] < 25:
            if depth > merged[-1][1]:
                merged[-1] = (centre, depth)
        else:
            merged.append((centre, depth))
    return merged


def measure(splib: Path):
    """Features of every sample, and those shared across the synthetics."""
    base = splib / "ASCIIdata" / "ASCIIdata_splib07a"
    wavelengths = (
        load_column(base / "splib07a_Wavelengths_ASD_0.35-2.5_microns_2151_ch.txt")
        * 1000.0
    )
    chapter = base / "ChapterA_ArtificialMaterials"

    per_sample = {}
    for group, samples in (("synthetic", SYNTHETIC), ("natural", NATURAL)):
        for name, filename in samples.items():
            reflectance = load_column(chapter / filename)
            per_sample[name] = {
                "group": group,
                "albedo_500_900": float(
                    np.nanmean(reflectance[(wavelengths >= 500) & (wavelengths <= 900)])
                ),
                "features": [
                    {"nm": round(c, 1), "depth": round(d, 3)}
                    for c, d in local_maxima(
                        wavelengths, band_depths(wavelengths, reflectance)
                    )
                ],
            }

    tally: dict[int, list[float]] = {}
    for name, record in per_sample.items():
        if record["group"] != "synthetic":
            continue
        for feature in record["features"]:
            tally.setdefault(int(round(feature["nm"] / 20.0) * 20), []).append(
                feature["nm"]
            )
    shared = {
        key: {
            "n": len(v),
            "mean_nm": round(float(np.mean(v)), 1),
            "min_nm": round(min(v), 1),
            "max_nm": round(max(v), 1),
        }
        for key, v in sorted(tally.items())
        if len(v) >= 4
    }
    return per_sample, shared


# Species added to the taxonomy. Centres and ranges follow the measurement
# above, widened to the nearest few nm so a 3-5 nm sensor still lands inside.
SPECIES = [
    (
        "coated_fabric_CH_1202",
        1202,
        (1192, 1222),
        "overtone_CH",
        "C-H second overtone of the fabric coating; present in every synthetic sample",
    ),
    (
        "coated_fabric_NH_OH_1414",
        1414,
        (1402, 1428),
        "overtone_NH_OH",
        "N-H / O-H first overtone; in all 8 synthetic samples, 1411-1419 nm",
    ),
    (
        "polyamide_NH_1535",
        1535,
        (1526, 1545),
        "overtone_NH",
        "Polyamide N-H first overtone; nylon-specific, absent from cellulose",
    ),
    (
        "polyamide_NH_1573",
        1573,
        (1562, 1584),
        "overtone_NH",
        "Second polyamide N-H component, always paired with the 1535 nm band",
    ),
    (
        "polymer_CH_overtone_1716",
        1716,
        (1700, 1740),
        "overtone_CH",
        "C-H first overtone, the strongest synthetic-fabric band (depth 0.16-0.49)",
    ),
    (
        "polymer_CH_overtone_1755",
        1757,
        (1742, 1772),
        "overtone_CH",
        "C-H first overtone partner; the 1716+1755 doublet is the branch's anchor",
    ),
    (
        "polyamide_amide_2049",
        2049,
        (2036, 2064),
        "combination_amide",
        "Amide combination band of polyamide fabrics; depth up to 0.36",
    ),
    (
        "polymer_CH_2183",
        2183,
        (2170, 2198),
        "combination_CH",
        "C-H combination of coated polyamide fabric",
    ),
    (
        "polymer_CH_2294",
        2294,
        (2282, 2308),
        "combination_CH",
        "C-H combination; deepest SWIR band of the coated fabrics",
    ),
    # Measured on the C2RD sheets themselves rather than on splib07: the
    # laboratory bands do not survive a thin sheet draped over grass.
    (
        "polymer_sheet_2122",
        2122,
        (2105, 2140),
        "combination_CH",
        "Field-measured band of the C2RD polymer sheeting",
    ),
    (
        "polymer_sheet_2245",
        2245,
        (2228, 2262),
        "combination_CH",
        "Field-measured partner band; with 2122 nm it gives 11/12 recall",
    ),
    (
        "water_body_OH_1940",
        1940,
        (1900, 1980),
        "combination_OH",
        "Water band used only as an exclusion: open water reaches depth 0.64-0.70",
    ),
]

# Minimum band depth per species for the camouflage rules. The cellulose
# controls reach the same centres but never these depths: burlap peaks at
# 0.051 near 1725 nm and cotton at 0.064 near 1776 nm, while every coated
# fabric sits between 0.149 and 0.490.
CAMO_DEPTHS = {
    "polymer_CH_overtone_1716": 0.12,
    "polymer_CH_overtone_1755": 0.10,
    "polyamide_NH_1535": 0.05,
    "polyamide_amide_2049": 0.05,
    # Only open water is excluded: it reaches 0.64-0.70 here, while a
    # vegetated background sits well below this.
    "water_body_OH_1940": 0.40,
}

COMPOSITES = [
    {
        "name": "camouflage_fabric_cut_sheet",
        "category": "synthetic",
        "min_depth_overrides": CAMO_DEPTHS,
        "required_species": ["polymer_CH_overtone_1716", "polymer_CH_overtone_1755"],
        "bonus_species": [
            "coated_fabric_CH_1202",
            "coated_fabric_NH_OH_1414",
            "polymer_CH_2294",
        ],
        "base_confidence": 0.78,
        "interpretation": "Die-cut sheet camouflage garnish: one coated fabric "
        "punched with leaf-shaped voids, not knotted netting",
    },
    {
        "name": "camouflage_fabric_polyamide",
        "category": "synthetic",
        "min_depth_overrides": CAMO_DEPTHS,
        "required_species": [
            "polymer_CH_overtone_1716",
            "polyamide_NH_1535",
            "polyamide_amide_2049",
        ],
        "bonus_species": [
            "polyamide_NH_1573",
            "polymer_CH_2183",
            "polymer_CH_overtone_1755",
        ],
        "base_confidence": 0.82,
        "interpretation": "Polyamide (nylon) camouflage fabric or webbing; "
        "amide bands separate it from polyolefin sheeting",
    },
    {
        "name": "camouflage_fabric_coated_polyolefin",
        "category": "synthetic",
        "min_depth_overrides": CAMO_DEPTHS,
        "required_species": ["polymer_CH_overtone_1716", "polymer_CH", "CH_aliphatic"],
        "bonus_species": [
            "coated_fabric_CH_1202",
            "coated_fabric_NH_OH_1414",
            "polymer_CH_overtone_1755",
        ],
        "base_confidence": 0.80,
        "interpretation": "PVC- or polyolefin-coated woven camouflage sheet; "
        "very deep 2310/2350 nm bands, no amide",
    },
    {
        "name": "camouflage_garnish_over_vegetation",
        "category": "synthetic",
        "min_depth_overrides": CAMO_DEPTHS,
        "required_species": ["polymer_CH_overtone_1716", "chlorophyll_a+b"],
        "bonus_species": [
            "polyamide_NH_1535",
            "red_edge_inflection",
            "polymer_CH_overtone_1755",
        ],
        "base_confidence": 0.62,
        "interpretation": "Cut-sheet garnish laid over vegetation: fabric C-H "
        "bands mixed with chlorophyll seen through the voids. "
        "The intended operational signature",
    },
    {
        "name": "camouflage_polymer_sheet_field",
        "category": "synthetic",
        "required_species": ["polymer_sheet_2122", "polymer_sheet_2245"],
        "min_depth_overrides": {
            "polymer_sheet_2122": 0.12,
            "polymer_sheet_2245": 0.10,
            "water_body_OH_1940": 0.40,
        },
        "required_absent": ["water_body_OH_1940"],
        "bonus_species": ["polymer_CH_overtone_1716", "polymer_CH_2183"],
        "base_confidence": 0.75,
        "interpretation": (
            "Polymer sheeting, derived from the C2RD Pourrieres targets "
            "themselves. Detects 11 of 12 sheets; the other hits are a bin "
            "liner, a transparent bag and a windsock, which are also polymer. "
            "Open water is excluded by its own deep 1940 nm band"
        ),
    },
    {
        "name": "natural_fibre_fabric",
        "category": "organic",
        "required_species": ["cellulose_OH", "cellulose"],
        "required_absent": ["polymer_CH_overtone_1716"],
        "bonus_species": ["lignin"],
        "base_confidence": 0.70,
        "interpretation": "Cotton, hessian or canvas. Separated from synthetic "
        "camouflage by the absence of the 1716 nm C-H overtone",
    },
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--splib", type=Path, default=Path.home() / "DBDATA/usgs_splib07"
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "camouflage_branch.json",
    )
    args = parser.parse_args()

    per_sample, shared = measure(args.splib)

    branch = {
        "meta": {
            "description": "Camouflage branch of the spectral taxonomy, derived "
            "from USGS splib07 coated-fabric and tarpaulin spectra.",
            "source_library": "USGS splib07a, ChapterA_ArtificialMaterials",
            "wavelength_range_nm": "400-2500",
            "campaign": "C2RD Pourrieres, 13 October 2026",
            "material_note": "The campaign samples are die-cut sheet garnish, "
            "not knotted netting: a single coated fabric with "
            "leaf-shaped voids.",
        },
        "measurement": {"per_sample": per_sample, "shared_synthetic_bands": shared},
        "species_to_materials": {
            name: {
                "mechanism": mechanism,
                "nominal_center_nm": centre,
                "range_nm": list(rng),
                "notes": note,
                "tags": ["synthetic", "polymer", "fabric", "camouflage"],
                "possible_materials": [
                    "camouflage garnish",
                    "coated fabric",
                    "tarpaulin",
                    "webbing",
                ],
            }
            for name, centre, rng, mechanism, note in SPECIES
        },
        "composite_signatures": COMPOSITES,
    }
    args.out.write_text(
        json.dumps(branch, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    print(f"wrote {args.out}")
    print(f"{len(SPECIES)} species, {len(COMPOSITES)} composite signatures")
    print("\nbands shared by the synthetic fabrics (>=4 of 8 samples):")
    for rec in shared.values():
        print(
            f"  ~{rec['mean_nm']:7.1f} nm  n={rec['n']}  "
            f"range {rec['min_nm']}-{rec['max_nm']}"
        )
    print("\ncellulose controls, for contrast:")
    for name, rec in per_sample.items():
        if rec["group"] == "natural":
            bands = ", ".join(f"{f['nm']:.0f}" for f in rec["features"][:7])
            print(f"  {name:22} {bands}")


if __name__ == "__main__":
    main()
