#!/usr/bin/env python
##############################################################################
# MODULE:    i.hyper.spectroscopy
# AUTHOR(S): Spectral Feature Extraction and Interpretation Engine
# PURPOSE:   Per-pixel spectral interpretation of hyperspectral 3D raster:
#            assigns material/class hypotheses from absorption features
#            (400–2500 nm VNIR+SWIR) using physics-based spectral rules.
# COPYRIGHT: (C) 2025 by the GRASS Development Team
# SPDX-License-Identifier: GPL-2.0-or-later
##############################################################################

# %module
# % description: Per-pixel spectral feature interpretation of hyperspectral imagery (400–2500 nm)
# % keyword: imagery
# % keyword: hyperspectral
# % keyword: spectroscopy
# % keyword: classification
# % keyword: minerals
# % keyword: vegetation
# %end

# %option G_OPT_R3_INPUT
# % key: input
# % required: yes
# % description: Input hyperspectral 3D raster map
# % guisection: Input
# %end

# %option G_OPT_R_OUTPUT
# % key: output
# % required: no
# % description: Output dominant-class raster (integer class codes)
# % guisection: Output
# %end

# %option
# % key: confidence
# % type: string
# % required: no
# % description: Output confidence raster name (0–1)
# % guisection: Output
# %end

# %option
# % key: min_band_depth
# % type: double
# % required: no
# % answer: 0.02
# % description: Minimum band depth to count as detected absorption feature (0–1)
# % guisection: Processing
# %end

# %option
# % key: coordinates
# % type: string
# % required: no
# % description: East,North coordinates for single-point interpretation (point mode)
# % guisection: Point mode
# %end

# %flag
# % key: n
# % description: Only use bands marked as valid (valid=1) in metadata
# % guisection: Processing
# %end

# %flag
# % key: i
# % description: Print band/wavelength info and exit without processing
# % guisection: Processing
# %end

# %flag
# % key: p
# % description: Point mode: interpret spectrum at coordinates (requires coordinates=)
# % guisection: Point mode
# %end

# %flag
# % key: v
# % description: Verbose: print full hypothesis ranking for point mode
# % guisection: Point mode
# %end

from __future__ import annotations

import sys
import os
import ctypes
import atexit
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import grass.script as gs

# ---------------------------------------------------------------------------
# Temporary raster cleanup
# ---------------------------------------------------------------------------

_TMP_RASTERS: list[str] = []


def _cleanup_tmp_rasters():
    if _TMP_RASTERS:
        gs.run_command('g.remove', type='raster',
                       name=','.join(_TMP_RASTERS), flags='f', quiet=True)


atexit.register(_cleanup_tmp_rasters)

# ---------------------------------------------------------------------------
# Fast 3D raster Z-slice extraction (tile-bulk, RASTER3D_NO_CACHE)
# ---------------------------------------------------------------------------

_G3D_LIB = None


def _load_g3d_lib():
    """Load libgrass_g3d and set up Rast3d_extract_z_slice() signature."""
    grass_config = gs.read_command('grass', '--config', 'path').strip()
    lib_path = os.path.join(grass_config, 'lib', 'libgrass_g3d.so')
    lib = ctypes.CDLL(lib_path)
    lib.Rast3d_extract_z_slice.restype = ctypes.c_int
    lib.Rast3d_extract_z_slice.argtypes = [
        ctypes.c_char_p,   # name3d
        ctypes.c_char_p,   # mapset3d (NULL = search)
        ctypes.c_int,      # z  (0-based)
        ctypes.c_char_p,   # name2d
    ]
    return lib


def extract_band(raster3d: str, band_num: int) -> str:
    """Extract band_num (1-based) from raster3d into a temporary 2D raster.

    Uses Rast3d_extract_z_slice() with RASTER3D_NO_CACHE: each tile at the
    target Z level is read exactly once (tile-bulk path), vs. one function
    call per voxel with the default cache API.

    Returns the temporary 2D raster name (registered for cleanup at exit).
    """
    global _G3D_LIB
    if _G3D_LIB is None:
        _G3D_LIB = _load_g3d_lib()

    z = band_num - 1  # 1-based → 0-based
    base = raster3d.replace('@', '_').replace('#', '_').replace('.', '_')
    tmp_name = f"tmp_spectroscopy_{base}_{band_num}"

    if '@' in raster3d:
        name3d, mapset3d = raster3d.split('@', 1)
    else:
        name3d, mapset3d = raster3d, ''

    ret = _G3D_LIB.Rast3d_extract_z_slice(
        name3d.encode(),
        mapset3d.encode() if mapset3d else None,
        ctypes.c_int(z),
        tmp_name.encode(),
    )
    if ret != 0:
        gs.fatal(f"Rast3d_extract_z_slice failed for band {band_num} of {raster3d}")

    _TMP_RASTERS.append(tmp_name)
    return tmp_name

# ---------------------------------------------------------------------------
# Band metadata
# ---------------------------------------------------------------------------


def get_band_info(raster3d: str, only_valid: bool = False) -> list[dict]:
    """Return sorted list of {band, wavelength, fwhm, valid} dicts.

    Parses from r3.info -h history (format: 'Band N: W nm, FWHM: F nm').
    Falls back to r.support metadata per band if history is missing.
    """
    info = gs.raster3d_info(raster3d)
    depths = int(info['depths'])

    bands: list[dict] = []

    # Primary: parse r3.info history
    try:
        history = gs.read_command('r3.info', flags='h', map=raster3d)
        for line in history.split('\n'):
            line = line.strip()
            if not line.startswith('Band '):
                continue
            try:
                parts = line.split('Band ')[1].split(':')
                band_num = int(parts[0].strip())
                wavelength = float(parts[1].split('nm')[0].strip())
                fwhm_str = parts[2].split('nm')[0].strip() if len(parts) > 2 else '10'
                fwhm = float(fwhm_str) if fwhm_str else 10.0
                bands.append({'band': band_num, 'wavelength': wavelength,
                               'fwhm': fwhm, 'valid': True})
            except (ValueError, IndexError):
                pass
    except Exception:
        pass

    # Fallback: r.support per-band metadata
    if not bands:
        gs.verbose("No band info in r3.info history; trying r.support per band")
        for i in range(1, depths + 1):
            band_name = f"{raster3d}#{i}"
            wl = fwhm = None
            valid = True
            unit = 'nm'
            try:
                result = gs.read_command('r.support', map=band_name, flags='n')
                for ln in result.split('\n'):
                    ln = ln.strip()
                    if ln.startswith('wavelength='):
                        wl = float(ln.split('=')[1])
                    elif ln.startswith('FWHM='):
                        fwhm = float(ln.split('=')[1])
                    elif ln.startswith('valid='):
                        valid = int(ln.split('=')[1]) == 1
                    elif ln.startswith('unit='):
                        unit = ln.split('=')[1].strip()
            except Exception:
                pass
            if wl is None:
                continue
            # Unit conversion
            unit = unit.lower()
            if unit in ('um', 'µm', 'micrometer', 'micron'):
                wl *= 1000.0
            elif unit in ('m', 'meter'):
                wl *= 1e9
            bands.append({'band': i, 'wavelength': wl,
                          'fwhm': fwhm or 10.0, 'valid': valid})

    if not bands:
        gs.fatal(
            f"No wavelength metadata in '{raster3d}'. "
            "Import data with i.hyper.import or add wavelength metadata."
        )

    bands.sort(key=lambda b: b['wavelength'])

    if only_valid:
        bands = [b for b in bands if b['valid']]
        if not bands:
            gs.fatal("No valid bands found (all marked valid=0).")

    return bands

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class AbsorptionFeature:
    """A single absorption feature from the spectrum."""
    center_nm: float
    fwhm_nm: float
    reflectance: float          # reflectance at band bottom (0–1)
    depth: float = 0.0          # band depth vs local continuum (0–1)
    asymmetry: str = "unknown"  # "sharp" | "moderate" | "broad"
    feature_type: str = "absorption"


@dataclass
class ChromophoreAssignment:
    """A candidate absorbing species matching one or more features."""
    species: str
    mechanism: str
    supporting_features: list[float] = field(default_factory=list)
    confidence: float = 0.0
    notes: str = ""


@dataclass
class MaterialHypothesis:
    """A ranked candidate material / class."""
    name: str
    category: str              # mineral | biologic | organic | geological | synthetic
    confidence: float          # 0–1
    evidence: list[str] = field(default_factory=list)
    atomic_species: list[str] = field(default_factory=list)


@dataclass
class SpectralInterpretation:
    """Top-level result returned by interpret_spectrum()."""
    features: list[AbsorptionFeature]
    chromophores: list[ChromophoreAssignment]
    hypotheses: list[MaterialHypothesis]      # sorted descending by confidence
    dominant_class: Optional[str]
    atomic_summary: list[str]
    flag_notes: list[str]

# ---------------------------------------------------------------------------
# Spectral physics: absorption feature database
# ---------------------------------------------------------------------------

SPECTRAL_FEATURES_DB: list[dict] = [

    # ── ELECTRONIC / CHARGE-TRANSFER / CRYSTAL-FIELD ──────────────────────

    # Chlorophylls / pigments
    {"center": 430,  "range": (415, 445), "species": "chlorophyll_a+b",
     "mechanism": "π→π* Soret band", "tags": ["vegetation", "algae", "biologic"],
     "notes": "Blue Soret; paired with 670nm confirms chlorophyll"},
    {"center": 480,  "range": (460, 500), "species": "carotenoids",
     "mechanism": "π→π*", "tags": ["vegetation", "biologic", "pigment"],
     "notes": "β-carotene / xanthophyll blue-green absorption"},
    {"center": 550,  "range": (535, 565), "species": "chlorophyll_green_peak",
     "mechanism": "reflectance_peak", "tags": ["vegetation"],
     "notes": "Green reflectance peak; not an absorption"},
    {"center": 620,  "range": (605, 640), "species": "phycocyanin",
     "mechanism": "π→π*", "tags": ["cyanobacteria", "algae", "biologic"],
     "notes": "Dominant in cyanobacteria, absent in higher plants"},
    {"center": 665,  "range": (650, 680), "species": "chlorophyll_a",
     "mechanism": "π→π* Qy", "tags": ["vegetation", "algae", "biologic"],
     "notes": "Red Qy; depth correlates with LAI and chlorophyll content"},
    {"center": 680,  "range": (670, 695), "species": "chlorophyll_a+b",
     "mechanism": "π→π* Qy combined", "tags": ["vegetation", "biologic"],
     "notes": "Combined Chl-a+b; used in red-edge index"},
    {"center": 710,  "range": (700, 730), "species": "red_edge_inflection",
     "mechanism": "chlorophyll_edge", "tags": ["vegetation", "biologic"],
     "notes": "Inflection of red-edge; position shifts with chlorophyll stress"},
    {"center": 740,  "range": (730, 760), "species": "red_edge_plateau",
     "mechanism": "NIR_plateau_start", "tags": ["vegetation"],
     "notes": "Shoulder marking transition to NIR plateau"},

    # Hemoglobin / blood
    {"center": 415,  "range": (408, 425), "species": "oxyhemoglobin_soret",
     "mechanism": "π→π* Soret", "tags": ["blood", "biologic", "medical"],
     "notes": "Strong Soret band; differentiates blood from other pigments"},
    {"center": 542,  "range": (535, 550), "species": "oxyhemoglobin",
     "mechanism": "Q-band α", "tags": ["blood", "biologic"],
     "notes": "α Q-band HbO2; paired with 577nm"},
    {"center": 577,  "range": (570, 585), "species": "oxyhemoglobin",
     "mechanism": "Q-band β", "tags": ["blood", "biologic"],
     "notes": "β Q-band; doublet 542+577 = HbO2 signature"},
    {"center": 560,  "range": (550, 575), "species": "deoxyhemoglobin",
     "mechanism": "Q-band", "tags": ["blood", "biologic", "medical"],
     "notes": "Single broad band vs HbO2 doublet; indicates deoxygenation"},
    {"center": 630,  "range": (620, 640), "species": "methemoglobin",
     "mechanism": "Q-band", "tags": ["blood", "biologic", "medical"],
     "notes": "MetHb marker; clinically significant"},

    # Iron oxides / hydroxides
    {"center": 490,  "range": (460, 520), "species": "Fe3+_goethite",
     "mechanism": "crystal_field_6A1→4T2", "tags": ["mineral", "geological", "soil"],
     "notes": "Goethite spin-forbidden; intensity low but diagnostic"},
    {"center": 530,  "range": (510, 560), "species": "Fe3+_hematite",
     "mechanism": "crystal_field", "tags": ["mineral", "geological", "soil", "Mars"],
     "notes": "Hematite crystal-field; deep red color"},
    {"center": 700,  "range": (670, 740), "species": "Fe3+_charge_transfer",
     "mechanism": "ligand→metal CT", "tags": ["mineral", "geological"],
     "notes": "Broad shoulder; many Fe3+ oxides"},
    {"center": 900,  "range": (850, 960), "species": "Fe2+_olivine_pyroxene",
     "mechanism": "crystal_field_5T2→5E", "tags": ["mineral", "geological", "mafic"],
     "notes": "Band I; position depends on Mg/Fe ratio in olivine/pyroxene"},
    {"center": 1000, "range": (950, 1060), "species": "Fe2+_Band_I_orthopyroxene",
     "mechanism": "crystal_field M2", "tags": ["mineral", "geological", "mafic"],
     "notes": "OPx Band I; shifts with Ca content"},
    {"center": 1250, "range": (1190, 1310), "species": "Fe2+_olivine",
     "mechanism": "crystal_field M1+M2 combined", "tags": ["mineral", "geological"],
     "notes": "Olivine triplet shoulder region"},
    {"center": 2000, "range": (1950, 2060), "species": "Fe2+_Band_II_pyroxene",
     "mechanism": "crystal_field M2", "tags": ["mineral", "geological", "mafic"],
     "notes": "Band II; Band I/Band II ratio diagnoses pyroxene type"},

    # REE / rare earth
    {"center": 580,  "range": (575, 590), "species": "Nd3+",
     "mechanism": "4f intraconfigurational", "tags": ["REE", "mineral", "geological"],
     "notes": "Neodymium; narrow f→f transition diagnostic"},
    {"center": 745,  "range": (740, 755), "species": "Nd3+",
     "mechanism": "4f", "tags": ["REE", "mineral"],
     "notes": "Second Nd3+ band; paired with 580nm"},
    {"center": 800,  "range": (793, 812), "species": "Nd3+",
     "mechanism": "4f", "tags": ["REE", "mineral"],
     "notes": "Third Nd3+ feature"},
    {"center": 867,  "range": (860, 875), "species": "Nd3+",
     "mechanism": "4f", "tags": ["REE", "mineral"],
     "notes": "Fourth Nd3+ line"},
    {"center": 525,  "range": (518, 535), "species": "Sm3+_or_Er3+",
     "mechanism": "4f", "tags": ["REE", "mineral"],
     "notes": "Samarium or erbium; narrow"},
    {"center": 940,  "range": (930, 955), "species": "Pr3+",
     "mechanism": "4f", "tags": ["REE", "mineral"],
     "notes": "Praseodymium; rare"},

    # Copper / cobalt / chromium minerals
    {"center": 435,  "range": (420, 455), "species": "Co2+",
     "mechanism": "crystal_field", "tags": ["mineral", "pigment"],
     "notes": "Cobalt blue; smaltite / erythrite"},
    {"center": 700,  "range": (680, 730), "species": "Cu2+",
     "mechanism": "crystal_field Jahn-Teller", "tags": ["mineral", "geological"],
     "notes": "Turquoise, malachite, azurite broad Cu absorption"},
    {"center": 630,  "range": (610, 655), "species": "Cr3+",
     "mechanism": "crystal_field 4A2→4T2", "tags": ["mineral", "geological"],
     "notes": "Chromium; emerald, uvarovite, chrome-diopside"},
    {"center": 460,  "range": (440, 480), "species": "Cr3+",
     "mechanism": "crystal_field 4A2→4T1", "tags": ["mineral", "geological"],
     "notes": "Second Cr band; paired with 630"},
    {"center": 550,  "range": (535, 570), "species": "Mn2+_or_Mn3+",
     "mechanism": "crystal_field", "tags": ["mineral", "geological"],
     "notes": "Spessartine garnet, rhodochrosite, piemontite"},

    # Atmospheric O2 / H2O (reference / contamination check)
    {"center": 690,  "range": (685, 695), "species": "O2_A_band_atmospheric",
     "mechanism": "electronic_O2", "tags": ["atmospheric", "reference"],
     "notes": "Telluric O2 A-band; seen in solar-reflected spectra"},
    {"center": 760,  "range": (755, 765), "species": "O2_B_band_atmospheric",
     "mechanism": "electronic_O2", "tags": ["atmospheric", "reference"],
     "notes": "O2 B-band telluric"},
    {"center": 820,  "range": (810, 835), "species": "H2O_atmospheric",
     "mechanism": "rotational_vibrational", "tags": ["atmospheric", "reference"],
     "notes": "Water vapor column; band depth ~ pwv"},
    {"center": 940,  "range": (910, 970), "species": "H2O_atmospheric",
     "mechanism": "rotational_vibrational", "tags": ["atmospheric", "reference"],
     "notes": "Primary atmospheric water band"},
    {"center": 1140, "range": (1110, 1170), "species": "H2O_atmospheric",
     "mechanism": "rotational_vibrational", "tags": ["atmospheric", "reference"],
     "notes": "Secondary water vapor band"},

    # ── VIBRATIONAL OVERTONE / COMBINATION BANDS (SWIR) ───────────────────

    # Water (liquid / bound / ice)
    {"center": 970,  "range": (950, 990),  "species": "liquid_water_OH",
     "mechanism": "2nd_overtone_OH_stretch", "tags": ["water", "biologic", "mineral", "soil"],
     "notes": "2ν_OH; liquid water. Absent/shifted in ice"},
    {"center": 1200, "range": (1180, 1220), "species": "liquid_water_OH",
     "mechanism": "combination_bend+stretch", "tags": ["water", "biologic", "mineral"],
     "notes": "ν+δ combination; weaker than 1450"},
    {"center": 1450, "range": (1400, 1500), "species": "OH_water",
     "mechanism": "1st_overtone_OH_stretch", "tags": ["water", "mineral", "biologic", "soil"],
     "notes": "Strong; liquid water ~1450, structural OH shifts 1380–1420"},
    {"center": 1940, "range": (1900, 1980), "species": "liquid_water",
     "mechanism": "combination_stretch+bend", "tags": ["water", "biologic", "soil"],
     "notes": "H2O bend+stretch combination; liquid water diagnostic"},

    # Ice
    {"center": 1030, "range": (1000, 1060), "species": "ice_H2O",
     "mechanism": "OH_combination", "tags": ["ice", "cryosphere", "geological"],
     "notes": "Ice band shifts vs liquid water"},
    {"center": 1270, "range": (1240, 1290), "species": "ice_H2O",
     "mechanism": "OH_combination", "tags": ["ice", "cryosphere"],
     "notes": "Ice overtone; paired with 1030 confirms ice"},
    {"center": 1500, "range": (1480, 1530), "species": "ice_H2O",
     "mechanism": "1st_overtone_OH", "tags": ["ice", "cryosphere"],
     "notes": "Slightly red-shifted vs liquid water 1450"},
    {"center": 2000, "range": (1970, 2030), "species": "ice_H2O",
     "mechanism": "combination", "tags": ["ice", "cryosphere"],
     "notes": "Combination band; narrower than liquid"},

    # Hydroxyl minerals (clays, phyllosilicates)
    {"center": 1380, "range": (1360, 1410), "species": "structural_OH",
     "mechanism": "1st_overtone_OH_stretch", "tags": ["clay", "mineral", "geological"],
     "notes": "Bound OH overtone; kaolinite, serpentine, chlorite doublet ~1390+1410"},
    {"center": 1410, "range": (1395, 1430), "species": "Al-OH_clay",
     "mechanism": "OH_overtone", "tags": ["clay", "mineral", "geological"],
     "notes": "Al-OH; kaolinite, illite, muscovite"},
    {"center": 2200, "range": (2170, 2230), "species": "Al-OH",
     "mechanism": "combination_OH+Al-OH_bend", "tags": ["clay", "mineral", "geological", "altered"],
     "notes": "MOST DIAGNOSTIC clay feature; kaolinite 2200, smectite 2205, illite 2195–2210"},
    {"center": 2160, "range": (2140, 2175), "species": "NH4+",
     "mechanism": "combination_NH", "tags": ["mineral", "altered", "geological"],
     "notes": "Ammonium substitution in illite/muscovite"},
    {"center": 2250, "range": (2230, 2270), "species": "Mg-OH_Fe-OH",
     "mechanism": "combination", "tags": ["clay", "mineral", "geological"],
     "notes": "Mg-OH or Fe-OH; serpentine, chlorite, nontronite"},
    {"center": 2310, "range": (2290, 2340), "species": "Mg-OH_carbonate",
     "mechanism": "combination_Mg-OH+CO3", "tags": ["mineral", "geological", "carbonate"],
     "notes": "Dolomite, magnesite Mg-OH; also talc"},
    {"center": 2350, "range": (2330, 2380), "species": "carbonate_CO3",
     "mechanism": "combination_CO3", "tags": ["carbonate", "mineral", "geological"],
     "notes": "Calcite; paired with 2330 and 2500"},
    {"center": 2330, "range": (2310, 2355), "species": "carbonate_CO3",
     "mechanism": "combination_CO3_v1+v4", "tags": ["carbonate", "mineral", "geological"],
     "notes": "Calcite combination; doublet 2330+2350"},
    {"center": 2500, "range": (2470, 2510), "species": "carbonate_CO3",
     "mechanism": "combination_CO3_2nd", "tags": ["carbonate", "mineral", "geological"],
     "notes": "Third carbonate feature; confirms calcite / dolomite"},
    {"center": 2120, "range": (2095, 2145), "species": "carbonate_CO3",
     "mechanism": "combination", "tags": ["carbonate", "mineral", "geological"],
     "notes": "CO3 first combination; weaker than 2330"},

    # Sulphates
    {"center": 1450, "range": (1440, 1470), "species": "SO4_H2O_gypsum",
     "mechanism": "OH_combination_hydration", "tags": ["sulphate", "mineral", "geological", "evaporite"],
     "notes": "Gypsum structural water; doublet ~1445+1490"},
    {"center": 1490, "range": (1475, 1510), "species": "SO4_H2O_gypsum",
     "mechanism": "OH_combination_hydration", "tags": ["sulphate", "mineral", "geological"],
     "notes": "Gypsum doublet second component"},
    {"center": 1750, "range": (1730, 1770), "species": "gypsum_SO4",
     "mechanism": "SO4_combination", "tags": ["sulphate", "mineral", "geological"],
     "notes": "Gypsum SO4 band; rare outside evaporites"},
    {"center": 2210, "range": (2195, 2230), "species": "Al-OH_alunite",
     "mechanism": "Al-OH_combination", "tags": ["sulphate", "mineral", "altered"],
     "notes": "Alunite / jarosite; advanced argillic alteration"},
    {"center": 2265, "range": (2250, 2285), "species": "jarosite",
     "mechanism": "Fe-OH_combination", "tags": ["sulphate", "mineral", "altered"],
     "notes": "Jarosite Fe-OH; paired with Fe3+ VIS absorptions"},

    # Carbonates (extra)
    {"center": 1850, "range": (1830, 1880), "species": "carbonate_CO3",
     "mechanism": "overtone_CO3_v3", "tags": ["carbonate", "mineral", "geological"],
     "notes": "Weak CO3 overtone; confirms carbonate"},

    # Cellulose / lignin / dry vegetation
    {"center": 1480, "range": (1460, 1505), "species": "cellulose_OH",
     "mechanism": "OH_overtone_polysaccharide", "tags": ["vegetation", "dry_veg", "organic", "biologic"],
     "notes": "Cellulose/starch OH; dry plant material"},
    {"center": 1820, "range": (1800, 1840), "species": "cellulose_CH",
     "mechanism": "CH_combination", "tags": ["vegetation", "dry_veg", "organic"],
     "notes": "Cellulose CH combination; dry biomass indicator"},
    {"center": 2100, "range": (2080, 2120), "species": "starch",
     "mechanism": "CO+OH_combination", "tags": ["vegetation", "food", "organic", "biologic"],
     "notes": "Starch diagnostic; food quality / grain"},
    {"center": 2270, "range": (2255, 2290), "species": "cellulose",
     "mechanism": "CO+CC_combination", "tags": ["vegetation", "dry_veg", "organic"],
     "notes": "Cellulose; C-O-C ring mode"},
    {"center": 2340, "range": (2320, 2360), "species": "lignin",
     "mechanism": "CH_aromatic", "tags": ["vegetation", "dry_veg", "organic"],
     "notes": "Lignin aromatic CH; wood / bark"},

    # Proteins / nitrogen
    {"center": 2054, "range": (2040, 2070), "species": "protein_NH",
     "mechanism": "NH_combination", "tags": ["biologic", "protein", "food", "soil"],
     "notes": "Protein N-H; correlates with leaf N content"},
    {"center": 2170, "range": (2155, 2190), "species": "protein_NH",
     "mechanism": "NH_combination_2nd", "tags": ["biologic", "protein", "food"],
     "notes": "Secondary NH band; amino acid / peptide bond"},
    {"center": 1510, "range": (1495, 1530), "species": "protein_NH",
     "mechanism": "1st_overtone_NH", "tags": ["biologic", "protein"],
     "notes": "First NH overtone"},
    {"center": 1680, "range": (1660, 1700), "species": "protein_amide",
     "mechanism": "amide_II_combination", "tags": ["biologic", "protein"],
     "notes": "Amide II band combination; protein secondary structure"},

    # Lipids / waxes / hydrocarbons
    {"center": 1210, "range": (1190, 1230), "species": "CH_lipid",
     "mechanism": "2nd_overtone_CH_stretch", "tags": ["lipid", "organic", "biologic", "petroleum"],
     "notes": "CH2 2nd overtone; lipid / wax content"},
    {"center": 1390, "range": (1370, 1415), "species": "CH_aromatic",
     "mechanism": "2nd_overtone_CH_aromatic", "tags": ["organic", "petroleum", "coal"],
     "notes": "Aromatic CH overtone; PAH / coal"},
    {"center": 1720, "range": (1700, 1745), "species": "CH2_lipid",
     "mechanism": "1st_overtone_CH_stretch", "tags": ["lipid", "organic", "petroleum", "biologic"],
     "notes": "Fatty acid CH2; first overtone of C-H stretch ~3450 cm-1"},
    {"center": 1760, "range": (1740, 1780), "species": "CH3_lipid",
     "mechanism": "1st_overtone_CH3", "tags": ["lipid", "organic", "petroleum"],
     "notes": "Terminal CH3 overtone; petroleum / wax"},
    {"center": 2310, "range": (2290, 2330), "species": "CH2_CH3_hydrocarbon",
     "mechanism": "combination_CH", "tags": ["petroleum", "organic", "coal", "synthetic"],
     "notes": "Strong CH2/CH3 combination; crude oil, coal, plastics"},
    {"center": 2350, "range": (2330, 2375), "species": "CH_aliphatic",
     "mechanism": "combination", "tags": ["petroleum", "organic", "synthetic"],
     "notes": "Aliphatic CH; overlaps carbonate – use FWHM to distinguish"},

    # Soil organic matter / humus
    {"center": 1900, "range": (1875, 1925), "species": "soil_organic_matter",
     "mechanism": "OH+CO_combination", "tags": ["soil", "organic", "geological"],
     "notes": "SOM broad feature; correlates with organic carbon %"},

    # Plastics / polymers
    {"center": 1150, "range": (1130, 1175), "species": "polyethylene_CH",
     "mechanism": "combination_CH", "tags": ["synthetic", "plastic", "polymer"],
     "notes": "PE; HDPE/LDPE diagnostic"},
    {"center": 1215, "range": (1200, 1235), "species": "polypropylene_CH",
     "mechanism": "combination_CH3", "tags": ["synthetic", "plastic", "polymer"],
     "notes": "PP methyl group; differs from PE"},
    {"center": 1660, "range": (1645, 1675), "species": "nylon_NH",
     "mechanism": "NH_overtone", "tags": ["synthetic", "plastic", "polymer"],
     "notes": "Nylon amide NH; similar to protein but context differs"},
    {"center": 2310, "range": (2295, 2330), "species": "polymer_CH",
     "mechanism": "CH_combination", "tags": ["synthetic", "plastic", "polymer"],
     "notes": "CH combination; present in many polymers"},

    # Quartz / silicates
    {"center": 1380, "range": (1360, 1400), "species": "Si-OH_quartz",
     "mechanism": "OH_combination", "tags": ["mineral", "geological", "silicate"],
     "notes": "Silanol OH; quartz, opal, glass"},
    {"center": 2200, "range": (2185, 2220), "species": "Si-OH",
     "mechanism": "Si-OH_combination", "tags": ["mineral", "geological", "silicate"],
     "notes": "Si-OH; opal, amorphous silica, altered feldspar"},

    # Feldspar
    {"center": 1250, "range": (1220, 1285), "species": "Fe2+_feldspar",
     "mechanism": "crystal_field", "tags": ["mineral", "geological", "feldspar"],
     "notes": "Fe2+ in anorthite site; darkens plagioclase"},

    # Phlogopite / biotite (mica)
    {"center": 2350, "range": (2335, 2370), "species": "Mg-OH_phlogopite",
     "mechanism": "Mg-OH_combination", "tags": ["mineral", "geological", "mica"],
     "notes": "Phlogopite; Mg-rich mica"},
    {"center": 2230, "range": (2215, 2250), "species": "Fe-OH_biotite",
     "mechanism": "Fe-OH_combination", "tags": ["mineral", "geological", "mica"],
     "notes": "Biotite Fe-OH; darker mica"},

    # Amphiboles
    {"center": 2310, "range": (2295, 2330), "species": "Mg-OH_tremolite",
     "mechanism": "Mg-OH_combination", "tags": ["mineral", "geological", "amphibole"],
     "notes": "Tremolite / actinolite Mg-OH"},
    {"center": 2370, "range": (2355, 2390), "species": "Mg-OH_antigorite",
     "mechanism": "Mg-OH_combination", "tags": ["mineral", "geological", "serpentine"],
     "notes": "Antigorite / lizardite serpentine"},

    # Epidote / zoisite
    {"center": 2250, "range": (2235, 2270), "species": "Fe-OH_epidote",
     "mechanism": "Fe-OH_combination", "tags": ["mineral", "geological", "altered"],
     "notes": "Epidote; propylitic alteration indicator"},

    # Apatite / phosphate
    {"center": 970,  "range": (955, 985),  "species": "REE_apatite",
     "mechanism": "4f_REE_in_apatite", "tags": ["mineral", "geological", "phosphate"],
     "notes": "Nd in apatite; narrows above ambient background"},

    # Chlorite
    {"center": 2250, "range": (2235, 2270), "species": "Mg-OH_Fe-OH_chlorite",
     "mechanism": "Mg-OH+Fe-OH_combination", "tags": ["mineral", "geological", "clay"],
     "notes": "Chlorite; position shifts with Fe/Mg ratio"},
    {"center": 2340, "range": (2325, 2360), "species": "Mg-OH_chlorite",
     "mechanism": "Mg-OH_combination", "tags": ["mineral", "geological", "clay"],
     "notes": "Chlorite second feature"},

    # Zeolites
    {"center": 1440, "range": (1420, 1460), "species": "zeolite_H2O",
     "mechanism": "cage_H2O_OH", "tags": ["mineral", "geological", "zeolite"],
     "notes": "Channel water in zeolite; broad"},
    {"center": 1910, "range": (1895, 1930), "species": "zeolite_H2O",
     "mechanism": "H2O_combination", "tags": ["mineral", "geological", "zeolite"],
     "notes": "Zeolite H2O combination"},

    # Asbestos (chrysotile)
    {"center": 2320, "range": (2305, 2340), "species": "Mg-OH_chrysotile",
     "mechanism": "Mg-OH_combination", "tags": ["mineral", "geological", "hazard", "asbestos"],
     "notes": "Chrysotile serpentine; hazardous material flag"},

    # Coal / char / black carbon
    {"center": 670,  "range": (650, 690),  "species": "coal_aromatic",
     "mechanism": "π→π*", "tags": ["organic", "coal", "geological"],
     "notes": "Aromatic ring absorption; rank correlates with depth"},
    {"center": 1390, "range": (1375, 1410), "species": "coal_CH",
     "mechanism": "aromatic_CH_overtone", "tags": ["organic", "coal"],
     "notes": "Aromatic CH; rank indicator"},
    {"center": 1700, "range": (1680, 1730), "species": "coal_CH2",
     "mechanism": "1st_overtone_aliphatic_CH", "tags": ["organic", "coal"],
     "notes": "Aliphatic CH2 overtone; decreases with rank"},
    {"center": 2300, "range": (2285, 2320), "species": "coal_CH_combo",
     "mechanism": "CH_combination", "tags": ["organic", "coal"],
     "notes": "CH combination; diagnostic of coal"},

    # ── STRESS PIGMENTS & AQUATIC PHYTOPLANKTON ───────────────────────────

    # Anthocyanins / zeaxanthin (vegetation stress indicators)
    {"center": 540,  "range": (510, 575), "species": "anthocyanins",
     "mechanism": "π→π*_flavylium_cation", "tags": ["vegetation", "biologic", "stress"],
     "notes": "Leaf anthocyanin; stress/UV-protection/senescence pigment; broad; Gitelson et al. 2001 RSE"},
    {"center": 531,  "range": (525, 538), "species": "zeaxanthin_PRI",
     "mechanism": "xanthophyll_de-epoxidation", "tags": ["vegetation", "biologic", "stress"],
     "notes": "PRI 531nm; zeaxanthin increases under light/heat stress; very narrow vs carotenoid broadband; Gamon et al. 1992"},

    # Aquatic phytoplankton pigments
    {"center": 545,  "range": (530, 560), "species": "phycoerythrin",
     "mechanism": "π→π*_phycoerythrobilin", "tags": ["algae", "biologic", "aquatic"],
     "notes": "Phycoerythrin dominant band; red algae + PE-type cyanobacteria; doublet 495+545nm; Simis et al. 2005 RSE"},
    {"center": 495,  "range": (483, 510), "species": "phycoerythrin",
     "mechanism": "π→π*_phycourobilin", "tags": ["algae", "biologic", "aquatic"],
     "notes": "Phycoerythrin 2nd band (phycourobilin type); doublet 495+545nm = phycoerythrin signature"},
    {"center": 650,  "range": (640, 662), "species": "allophycocyanin",
     "mechanism": "π→π*_phycocyanobilin", "tags": ["cyanobacteria", "algae", "biologic", "aquatic"],
     "notes": "Allophycocyanin APC; 650nm distinct from phycocyanin 620nm; terminal energy-transfer pigment"},
    {"center": 490,  "range": (472, 508), "species": "fucoxanthin",
     "mechanism": "π→π*_carotenoid_diatom", "tags": ["algae", "biologic", "aquatic", "diatom"],
     "notes": "Fucoxanthin in diatoms and brown algae; paired with Chl-a 665nm confirms diatom bloom; Babin et al. 2003 JGR"},

    # ── IRON MINERALOGY EXPANSION ─────────────────────────────────────────

    # Ferrihydrite (poorly crystalline Fe3+ oxyhydroxide)
    {"center": 430,  "range": (400, 465), "species": "ferrihydrite",
     "mechanism": "crystal_field_Fe3+_LMCT", "tags": ["mineral", "geological", "soil", "mine_drainage"],
     "notes": "Ferrihydrite broad Fe3+ + LMCT; blueshifted vs hematite 530nm; AMD indicator; Sherman & Waite 1985"},
    {"center": 750,  "range": (700, 810), "species": "ferrihydrite",
     "mechanism": "crystal_field_Fe3+_2nd", "tags": ["mineral", "geological", "soil", "mine_drainage"],
     "notes": "Second ferrihydrite feature; 430+750nm pair distinguishes from goethite/hematite; Schwertmann & Cornell 2000"},

    # Lepidocrocite γ-FeOOH
    {"center": 475,  "range": (455, 498), "species": "lepidocrocite",
     "mechanism": "crystal_field_6A1→4T2", "tags": ["mineral", "geological", "soil"],
     "notes": "γ-FeOOH lepidocrocite; blueshifted ~15nm vs goethite 490nm; waterlogged/gleyed soils; Scheinost et al. 1998"},
    {"center": 1000, "range": (985, 1020), "species": "lepidocrocite",
     "mechanism": "crystal_field_Fe3+_NIR", "tags": ["mineral", "geological", "soil"],
     "notes": "Lepidocrocite NIR feature; 475+1000nm pair confirms lepidocrocite over goethite"},

    # Siderite FeCO3
    {"center": 1080, "range": (980, 1180), "species": "siderite_Fe2+",
     "mechanism": "crystal_field_5T2→5E_carbonate", "tags": ["mineral", "geological", "carbonate"],
     "notes": "Siderite FeCO3 very broad Fe2+ band FWHM ~200nm; reducing/anoxic sedimentary facies; Hunt & Salisbury 1971"},
    {"center": 2340, "range": (2323, 2358), "species": "siderite_CO3",
     "mechanism": "CO3_combination_Fe", "tags": ["mineral", "geological", "carbonate"],
     "notes": "Siderite CO3 combination; slightly shifted vs calcite; paired with broad Fe2+ confirms siderite"},

    # ── RARE EARTH ELEMENTS EXPANSION ────────────────────────────────────

    # Dysprosium Dy3+
    {"center": 910,  "range": (903, 920), "species": "Dy3+",
     "mechanism": "4f_intraconfigurational", "tags": ["REE", "mineral", "geological"],
     "notes": "Dy3+ HREE; narrow 4f at 910nm; verify narrow FWHM vs broad H2O atmospheric; Rowan et al. 1986 Econ. Geol."},
    {"center": 1260, "range": (1252, 1271), "species": "Dy3+",
     "mechanism": "4f", "tags": ["REE", "mineral"],
     "notes": "Second Dy3+ band; 910+1260nm doublet = dysprosium signature; Turner et al. 2014 Geophysics"},

    # Ytterbium Yb3+
    {"center": 975,  "range": (968, 983), "species": "Yb3+",
     "mechanism": "4f_2F7/2→2F5/2", "tags": ["REE", "mineral", "geological"],
     "notes": "Yb3+ very narrow 4f; FWHM <15nm distinguishes from broad H2O atmospheric at 970nm; Rowan et al. 1986"},

    # Erbium Er3+ (additional bands beyond Sm3+_or_Er3+ at 525nm)
    {"center": 650,  "range": (642, 660), "species": "Er3+",
     "mechanism": "4f", "tags": ["REE", "mineral"],
     "notes": "Er3+ VIS band 650nm; narrow FWHM <20nm vs chlorophyll_a broad Qy; Cloutis et al. 2002 JGR"},
    {"center": 975,  "range": (968, 983), "species": "Er3+",
     "mechanism": "4f", "tags": ["REE", "mineral"],
     "notes": "Er3+ NIR; paired with 650nm and 1530nm confirms erbium over Yb3+"},
    {"center": 1530, "range": (1520, 1545), "species": "Er3+",
     "mechanism": "4f_4I13/2→4I15/2", "tags": ["REE", "mineral"],
     "notes": "Er3+ 1530nm telecom band; highly diagnostic; unique in SWIR; Cloutis et al. 2002"},

    # Holmium Ho3+
    {"center": 450,  "range": (442, 460), "species": "Ho3+",
     "mechanism": "4f_intraconfigurational", "tags": ["REE", "mineral"],
     "notes": "Ho3+ first 4f band; doublet with 537nm; USGS Spectral Library v7"},
    {"center": 537,  "range": (529, 546), "species": "Ho3+",
     "mechanism": "4f", "tags": ["REE", "mineral"],
     "notes": "Ho3+ second band; 450+537nm doublet = holmium signature"},

    # Europium Eu3+
    {"center": 535,  "range": (528, 543), "species": "Eu3+",
     "mechanism": "4f_7F0→5D1", "tags": ["REE", "mineral"],
     "notes": "Eu3+ very narrow 8nm band; requires high-resolution sensor (>5nm sampling); USGS Spec. Lib. v7"},

    # Terbium Tb3+
    {"center": 490,  "range": (483, 498), "species": "Tb3+",
     "mechanism": "4f_intraconfigurational", "tags": ["REE", "mineral"],
     "notes": "Tb3+ narrow 4f; overlaps Fe3+ goethite and fucoxanthin — requires FWHM <12nm; Turner et al. 2014"},

    # Cerium Ce3+
    {"center": 400,  "range": (390, 435), "species": "Ce3+",
     "mechanism": "4f→5d_charge_transfer", "tags": ["REE", "mineral", "geological"],
     "notes": "Ce3+ 4f→5d CT; broad UV edge; band onset in 400–430nm; Ce-rich carbonatite/monazite; Rowan et al. 1986"},

    # ── CLAY / PHYLLOSILICATE EXPANSION ──────────────────────────────────

    # Pyrophyllite Al2Si4O10(OH)2
    {"center": 2165, "range": (2150, 2181), "species": "Al-OH_pyrophyllite",
     "mechanism": "Al-OH_combination", "tags": ["clay", "mineral", "geological", "altered"],
     "notes": "Pyrophyllite 2165nm; 35nm BLUESHIFTED vs kaolinite 2200nm; epithermal acid-sulfate; Chabrillat et al. 2002"},
    {"center": 1390, "range": (1378, 1403), "species": "Al-OH_pyrophyllite",
     "mechanism": "Al-OH_overtone", "tags": ["clay", "mineral", "geological"],
     "notes": "Pyrophyllite 1390nm OH overtone; narrower and blueshifted vs kaolinite Al-OH_clay 1410nm"},

    # Gibbsite Al(OH)3
    {"center": 2263, "range": (2248, 2279), "species": "gibbsite",
     "mechanism": "Al-OH_combination", "tags": ["mineral", "geological", "soil", "laterite"],
     "notes": "Gibbsite Al(OH)3 2263nm; distinct from kaolinite 2200nm; tropical laterite indicator; Clark et al. 1990"},
    {"center": 2387, "range": (2372, 2403), "species": "gibbsite",
     "mechanism": "Al-OH_combination_2nd", "tags": ["mineral", "geological", "soil", "laterite"],
     "notes": "Second gibbsite band; 2263+2387nm doublet = gibbsite signature; deep weathering profile mapping"},

    # Boehmite / diaspore AlOOH
    {"center": 2175, "range": (2160, 2191), "species": "boehmite",
     "mechanism": "Al-OH_combination", "tags": ["mineral", "geological", "bauxite"],
     "notes": "Boehmite/diaspore AlOOH 2175nm; distinct from gibbsite 2263nm and kaolinite 2200nm; USGS Spec. Lib. v7"},
    {"center": 2110, "range": (2095, 2128), "species": "boehmite",
     "mechanism": "Al-OH_combination_2nd", "tags": ["mineral", "geological", "bauxite"],
     "notes": "Second boehmite feature; doublet 2175+2110nm confirms boehmite over kaolinite; bauxite horizon"},

    # Halloysite Al2Si2O5(OH)4 (hydrated kaolinite)
    {"center": 2205, "range": (2190, 2221), "species": "halloysite",
     "mechanism": "Al-OH_combination", "tags": ["clay", "mineral", "geological"],
     "notes": "Halloysite Al-OH 2205nm; slightly shifted/broader vs kaolinite 2200nm; interlayer water diagnostic"},

    # Phengite (Si-rich muscovite)
    {"center": 2215, "range": (2205, 2228), "species": "phengite_Al-OH",
     "mechanism": "Al-OH_combination_Si-rich", "tags": ["mineral", "geological", "mica"],
     "notes": "Phengite 2215nm; shifted vs muscovite 2200nm by Si→Al substitution; porphyry ore pathfinder; Duke 1994"},

    # Palygorskite / attapulgite
    {"center": 2200, "range": (2183, 2218), "species": "palygorskite",
     "mechanism": "Al-Mg-OH_combination", "tags": ["clay", "mineral", "geological"],
     "notes": "Palygorskite/attapulgite Al-Mg-OH fibrous clay; arid climate paleosol indicator; Chabrillat et al. 2002"},
    {"center": 1415, "range": (1403, 1432), "species": "palygorskite",
     "mechanism": "structural_H2O_channel", "tags": ["clay", "mineral", "geological"],
     "notes": "Palygorskite channel water 1415nm; slightly shifted from kaolinite Al-OH_clay 1410nm"},

    # Sepiolite Mg4Si6O15(OH)2
    {"center": 2310, "range": (2294, 2328), "species": "sepiolite_Mg-OH",
     "mechanism": "Mg-OH_combination", "tags": ["clay", "mineral", "geological"],
     "notes": "Sepiolite Mg-OH fibrous clay 2310nm; arid lacustrine; Chabrillat et al. 2002"},
    {"center": 1450, "range": (1433, 1470), "species": "sepiolite_H2O",
     "mechanism": "zeolitic_H2O_adsorbed", "tags": ["clay", "mineral", "geological"],
     "notes": "Sepiolite zeolitic/adsorbed water 1450nm; broader than gypsum doublet"},

    # Talc Mg3Si4O10(OH)2
    {"center": 2315, "range": (2298, 2333), "species": "talc_Mg-OH",
     "mechanism": "Mg-OH_combination", "tags": ["mineral", "geological", "phyllosilicate"],
     "notes": "Talc Mg-OH 2315nm; 2315+2390nm doublet unique to talc; ultramafic talc-carbonate zones; Clark et al. 1990"},
    {"center": 2390, "range": (2375, 2407), "species": "talc_Mg-OH",
     "mechanism": "Mg-OH_combination_2nd", "tags": ["mineral", "geological", "phyllosilicate"],
     "notes": "Second talc band 2390nm; doublet 2315+2390nm is diagnostic of talc"},
    {"center": 1390, "range": (1373, 1407), "species": "talc_OH",
     "mechanism": "Mg-OH_overtone", "tags": ["mineral", "geological", "phyllosilicate"],
     "notes": "Talc Mg-OH overtone 1390nm; bonus confirmatory band"},

    # ── SULPHATE MINERAL EXPANSION ────────────────────────────────────────

    # Epsomite MgSO4·7H2O
    {"center": 1450, "range": (1432, 1473), "species": "epsomite_H2O",
     "mechanism": "H2O_overtone_structural", "tags": ["sulphate", "mineral", "geological", "evaporite"],
     "notes": "Epsomite MgSO4·7H2O structural water 1450nm; broader than gypsum; Mg-SO4 band required; Cloutis et al. 2006"},
    {"center": 1940, "range": (1918, 1963), "species": "epsomite_H2O",
     "mechanism": "H2O_combination_structural", "tags": ["sulphate", "mineral", "geological"],
     "notes": "Second epsomite water band 1940nm; paired with 1450nm"},
    {"center": 2200, "range": (2180, 2222), "species": "epsomite_MgSO4",
     "mechanism": "SO4_combination_Mg", "tags": ["sulphate", "mineral", "geological"],
     "notes": "Epsomite Mg-SO4 combination ~2200nm; water bands distinguish from Al-OH clays"},

    # Bassanite / anhydrite CaSO4·½H2O
    {"center": 2165, "range": (2148, 2183), "species": "bassanite_SO4",
     "mechanism": "SO4_combination_Ca", "tags": ["sulphate", "mineral", "geological", "evaporite"],
     "notes": "Bassanite/anhydrite 2165nm; absent strong H2O bands distinguishes from gypsum; Bishop & Murad 2004"},

    # Kieserite MgSO4·H2O
    {"center": 1640, "range": (1624, 1658), "species": "kieserite_H2O",
     "mechanism": "H2O_combination_monohydrate", "tags": ["sulphate", "mineral", "geological"],
     "notes": "Kieserite monohydrate 1640nm; distinct from epsomite heptahydrate 1450nm; Cloutis et al. 2006"},
    {"center": 2220, "range": (2205, 2238), "species": "kieserite_MgSO4",
     "mechanism": "SO4_combination_Mg", "tags": ["sulphate", "mineral", "geological"],
     "notes": "Kieserite Mg-SO4 2220nm; paired with 1640nm H2O confirms kieserite"},

    # Szomolnokite FeSO4·H2O
    {"center": 1000, "range": (978, 1025), "species": "szomolnokite_Fe2+",
     "mechanism": "crystal_field_Fe2+_SO4", "tags": ["sulphate", "mineral", "geological", "mine_drainage"],
     "notes": "Szomolnokite FeSO4·H2O Fe2+ ~1000nm; AMD iron sulfate mineral; Cloutis et al. 2006"},
    {"center": 1490, "range": (1472, 1512), "species": "szomolnokite_H2O",
     "mechanism": "H2O_combination_FeSO4", "tags": ["sulphate", "mineral", "geological"],
     "notes": "Szomolnokite structural H2O 1490nm; paired with Fe2+ 1000nm confirms szomolnokite"},

    # ── ADDITIONAL GEOLOGICAL MINERALS ───────────────────────────────────

    # Prehnite Ca2Al2Si3O10(OH)2
    {"center": 2243, "range": (2230, 2258), "species": "prehnite_Al-OH",
     "mechanism": "Al-OH_Ca-OH_combination", "tags": ["mineral", "geological", "altered", "metamorphic"],
     "notes": "Prehnite 2243nm; diagnostic vs epidote 2250nm and chlorite 2250nm; low-grade metamorphic; Kruse et al. 2012 RSE"},
    {"center": 1490, "range": (1477, 1508), "species": "prehnite_OH",
     "mechanism": "OH_combination", "tags": ["mineral", "geological"],
     "notes": "Prehnite structural OH 1490nm; paired with 2243nm confirms prehnite"},

    # Tourmaline / schorl (boron silicate)
    {"center": 2120, "range": (2103, 2140), "species": "schorl_B-OH",
     "mechanism": "B-OH_combination_tourmaline", "tags": ["mineral", "geological", "tourmaline"],
     "notes": "Schorl/tourmaline B-OH combination 2120nm; unique to boron silicates; metasediment/greisen indicator"},

    # Crocidolite / riebeckite (blue asbestos)
    {"center": 2320, "range": (2302, 2342), "species": "crocidolite_Fe-OH",
     "mechanism": "Fe-OH_combination_amphibole", "tags": ["mineral", "geological", "hazard", "asbestos", "amphibole"],
     "notes": "Crocidolite riebeckite Fe-OH 2320nm; HAZARD blue asbestos; distinct from chrysotile by Fe3+ VIS band"},
    {"center": 440,  "range": (420, 463), "species": "crocidolite_Fe3+",
     "mechanism": "crystal_field_Fe3+_Na-amphibole", "tags": ["mineral", "geological", "hazard", "asbestos"],
     "notes": "Crocidolite Fe3+ in Na-amphibole 440nm; broad; paired with Fe-OH 2320nm confirms crocidolite"},

    # Ankerite Ca(Mg,Fe)(CO3)2
    {"center": 2320, "range": (2303, 2339), "species": "ankerite_Mg-OH",
     "mechanism": "Mg-OH_combination_carbonate", "tags": ["mineral", "geological", "carbonate"],
     "notes": "Ankerite Mg-OH 2320nm; intermediate between dolomite 2310nm and calcite; ore alteration halo; Clark et al. 1990"},
    {"center": 2335, "range": (2318, 2352), "species": "ankerite_CO3",
     "mechanism": "CO3_combination_ankerite", "tags": ["mineral", "geological", "carbonate"],
     "notes": "Ankerite CO3 combination 2335nm; doublet 2320+2335nm diagnostic of ankerite"},

    # ── PLASTIC POLYMER EXPANSION ─────────────────────────────────────────

    # PET (polyethylene terephthalate)
    {"center": 1730, "range": (1713, 1750), "species": "PET_CH",
     "mechanism": "CH_1st_overtone_aromatic+aliphatic", "tags": ["synthetic", "plastic", "polymer"],
     "notes": "PET 1730nm; aromatic+aliphatic CH; slightly redshifted vs PE 1720nm; Garaba & Dierssen 2018 Sci. Rep."},
    {"center": 2170, "range": (2155, 2188), "species": "PET_aromatic",
     "mechanism": "aromatic_CH_combination", "tags": ["synthetic", "plastic", "polymer"],
     "notes": "PET aromatic CH combination 2170nm; paired with 1730nm confirms PET over PE/PP"},

    # PVC (polyvinyl chloride)
    {"center": 1730, "range": (1712, 1750), "species": "PVC_CH2",
     "mechanism": "CH2_1st_overtone_PVC", "tags": ["synthetic", "plastic", "polymer"],
     "notes": "PVC CH2 overtone 1730nm; distinguish from PET by unique C-Cl band at 2030nm; Garaba & Dierssen 2018"},
    {"center": 2030, "range": (2010, 2058), "species": "PVC_CCl",
     "mechanism": "C-Cl_combination", "tags": ["synthetic", "plastic", "polymer"],
     "notes": "PVC C-Cl combination 2030nm; unique to chlorinated plastic; most diagnostic PVC feature"},

    # Polystyrene
    {"center": 1680, "range": (1667, 1695), "species": "polystyrene_CH_aromatic",
     "mechanism": "aromatic_CH_1st_overtone", "tags": ["synthetic", "plastic", "polymer"],
     "notes": "Polystyrene aromatic C-H first overtone 1680nm; narrower than protein_amide II; Li et al. 2022 RS"},
    {"center": 2240, "range": (2227, 2255), "species": "polystyrene_aromatic",
     "mechanism": "aromatic_ring_combination", "tags": ["synthetic", "plastic", "polymer"],
     "notes": "Polystyrene aromatic ring combination 2240nm; diagnostic narrow band; paired with 1680nm"},

    # Rubber / polyisoprene
    {"center": 1640, "range": (1622, 1658), "species": "rubber_C=C",
     "mechanism": "C=C_combination_isoprene", "tags": ["synthetic", "organic", "rubber"],
     "notes": "Polyisoprene C=C combination 1640nm; distinguishes rubber from saturated hydrocarbons"},
    {"center": 1720, "range": (1703, 1740), "species": "rubber_CH",
     "mechanism": "CH2_1st_overtone_rubber", "tags": ["synthetic", "organic", "rubber"],
     "notes": "Rubber CH2 overtone 1720nm; paired with C=C 1640nm confirms natural/synthetic rubber"},

    # ── URBAN MATERIALS ───────────────────────────────────────────────────

    # Concrete / hardened cement (portlandite Ca(OH)2)
    {"center": 1460, "range": (1443, 1479), "species": "portlandite_CaOH",
     "mechanism": "Ca-OH_overtone", "tags": ["urban", "mineral", "concrete"],
     "notes": "Portlandite Ca(OH)2 in hardened cement; 1460nm Ca-OH; paired with CO3 carbonation band; Herold et al. 2004"},

    # ── SOIL NUTRIENTS ────────────────────────────────────────────────────

    # Nitrate NO3-
    {"center": 2050, "range": (2030, 2075), "species": "nitrate_NO3",
     "mechanism": "NO3_combination_overtone", "tags": ["soil", "agricultural", "geochemical"],
     "notes": "Soil nitrate NO3- combination 2050nm; precision agriculture / geochemical survey; Malley et al. 1999 SSSAJ"},
    {"center": 1410, "range": (1397, 1425), "species": "nitrate_NO3",
     "mechanism": "NO3_combination", "tags": ["soil", "agricultural"],
     "notes": "Second nitrate band 1410nm; weaker; overlaps Al-OH_clay — paired 2050nm band required"},

    # Ammonium NH4+ second band (extends primary NH4+ at 2160nm)
    {"center": 1560, "range": (1546, 1577), "species": "NH4+_2nd",
     "mechanism": "NH4_combination_2nd", "tags": ["mineral", "soil", "altered"],
     "notes": "Second ammonium band 1560nm; paired with NH4+ 2160nm confirms ammonium-illite/buddingtonite; Krohn & Altaner 1987"},
]

# ---------------------------------------------------------------------------
# Composite rules: multi-species signatures → named material classes
# ---------------------------------------------------------------------------

COMPOSITE_RULES: list[dict] = [
    {"name": "green_vegetation",       "category": "biologic",
     "required_species": ["chlorophyll_a+b", "red_edge_inflection"],
     "bonus_species": ["carotenoids", "liquid_water_OH"], "base_score": 0.90,
     "notes": "Chlorophyll + red-edge = live green vegetation"},
    {"name": "stressed_vegetation",    "category": "biologic",
     "required_species": ["chlorophyll_a+b"], "required_absent": ["red_edge_inflection"],
     "bonus_species": ["cellulose_OH", "liquid_water_OH"], "base_score": 0.70,
     "notes": "Chl without strong red-edge: stressed or senescing plant"},
    {"name": "dry_plant_litter",       "category": "biologic",
     "required_species": ["cellulose_OH", "cellulose"],
     "bonus_species": ["lignin", "CH_lipid"], "base_score": 0.80,
     "notes": "Cellulose+lignin without liquid water: senescent/dry biomass"},
    {"name": "oxyhemoglobin_blood",    "category": "biologic",
     "required_species": ["oxyhemoglobin", "oxyhemoglobin_soret"],
     "base_score": 0.92, "notes": "Soret + Q-band doublet: oxygenated blood"},
    {"name": "deoxyhemoglobin_blood",  "category": "biologic",
     "required_species": ["deoxyhemoglobin", "oxyhemoglobin_soret"],
     "base_score": 0.85, "notes": "Soret + single Q: deoxygenated blood"},
    {"name": "kaolinite",              "category": "mineral",
     "required_species": ["Al-OH_clay", "Al-OH"],
     "base_score": 0.88, "notes": "1410+2200 doublet = kaolinite / halloysite"},
    {"name": "smectite_montmorillonite","category": "mineral",
     "required_species": ["Al-OH", "liquid_water_OH"],
     "bonus_species": ["OH_water"], "base_score": 0.75,
     "notes": "Smectite: 2205 + high water content"},
    {"name": "muscovite_illite",        "category": "mineral",
     "required_species": ["Al-OH_clay", "Al-OH"],
     "bonus_species": ["structural_OH"], "base_score": 0.78,
     "notes": "Muscovite/illite: 1410 + 2195–2210 nm"},
    {"name": "calcite",                "category": "mineral",
     "required_species": ["carbonate_CO3", "carbonate_CO3"],
     "base_score": 0.85, "notes": "Carbonate doublet 2330+2350 nm ± 2500"},
    {"name": "dolomite",               "category": "mineral",
     "required_species": ["Mg-OH_carbonate", "carbonate_CO3"],
     "base_score": 0.83, "notes": "Mg-OH at 2310 + carbonate features"},
    {"name": "gypsum_evaporite",       "category": "mineral",
     "required_species": ["SO4_H2O_gypsum", "SO4_H2O_gypsum"],
     "bonus_species": ["gypsum_SO4"], "base_score": 0.87,
     "notes": "Gypsum: doublet ~1445+1490 + optional 1750"},
    {"name": "liquid_water",           "category": "material",
     "required_species": ["liquid_water", "liquid_water_OH"],
     "base_score": 0.92, "notes": "1450 + 1940 = liquid H2O"},
    {"name": "ice_snow",               "category": "geological",
     "required_species": ["ice_H2O", "ice_H2O"],
     "base_score": 0.88, "notes": "1030 + 1270 + 1500 combination = H2O ice"},
    {"name": "olivine",                "category": "mineral",
     "required_species": ["Fe2+_olivine_pyroxene", "Fe2+_olivine"],
     "base_score": 0.80, "notes": "Broad Band I ~900 + shoulder 1250: olivine"},
    {"name": "pyroxene_orthopyroxene", "category": "mineral",
     "required_species": ["Fe2+_Band_I_orthopyroxene", "Fe2+_Band_II_pyroxene"],
     "base_score": 0.82, "notes": "Dual-band pyroxene; OPx Band I + Band II"},
    {"name": "hematite_iron_oxide",    "category": "mineral",
     "required_species": ["Fe3+_hematite", "Fe3+_charge_transfer"],
     "base_score": 0.84, "notes": "Hematite: VIS crystal-field + broad CT edge"},
    {"name": "goethite",               "category": "mineral",
     "required_species": ["Fe3+_goethite", "Fe3+_charge_transfer"],
     "base_score": 0.80, "notes": "Goethite: ~490 + broad Fe3+ CT"},
    {"name": "neodymium_REE",          "category": "mineral",
     "required_species": ["Nd3+", "Nd3+"],
     "base_score": 0.90, "notes": "Multiple Nd3+ 4f lines 580+745+800+867"},
    {"name": "crude_oil_petroleum",    "category": "organic",
     "required_species": ["CH_lipid", "CH2_lipid", "CH2_CH3_hydrocarbon"],
     "base_score": 0.86, "notes": "CH overtone series 1210+1720+2310"},
    {"name": "coal",                   "category": "geological",
     "required_species": ["coal_aromatic", "coal_CH_combo"],
     "bonus_species": ["coal_CH2"], "base_score": 0.82,
     "notes": "Aromatic VIS + SWIR CH combination"},
    {"name": "polyethylene_plastic",   "category": "synthetic",
     "required_species": ["polyethylene_CH", "polymer_CH"],
     "base_score": 0.84, "notes": "PE: 1150+2310 nm"},
    {"name": "protein_rich_material",  "category": "biologic",
     "required_species": ["protein_NH", "protein_NH"],
     "bonus_species": ["protein_amide"], "base_score": 0.78,
     "notes": "NH combination features: protein content"},
    {"name": "cyanobacteria_algae",    "category": "biologic",
     "required_species": ["phycocyanin", "chlorophyll_a+b"],
     "base_score": 0.87, "notes": "Phycocyanin + Chl: cyanobacteria or mixed algae"},
    {"name": "serpentinite",           "category": "mineral",
     "required_species": ["Mg-OH_antigorite", "Mg-OH_Fe-OH"],
     "base_score": 0.82, "notes": "Serpentine group: antigorite / lizardite"},
    {"name": "advanced_argillic_alteration", "category": "geological",
     "required_species": ["Al-OH_alunite", "Al-OH"],
     "base_score": 0.83, "notes": "Alunite/kaolinite: hydrothermal advanced argillic"},
    {"name": "propylitic_alteration",  "category": "geological",
     "required_species": ["Fe-OH_epidote", "Mg-OH_Fe-OH_chlorite"],
     "base_score": 0.80, "notes": "Epidote + chlorite: propylitic hydrothermal zone"},
    {"name": "chrysotile_asbestos",    "category": "mineral",
     "required_species": ["Mg-OH_chrysotile", "Mg-OH_antigorite"],
     "base_score": 0.85, "notes": "HAZARD: chrysotile asbestos fingerprint"},
    {"name": "starch_food_grain",      "category": "biologic",
     "required_species": ["starch", "cellulose_OH"],
     "base_score": 0.79, "notes": "Starch 2100 + cellulose 1480: grain / food"},

    # ── STRESS PIGMENTS & AQUATIC ────────────────────────────────────────
    {"name": "anthocyanin_stress",       "category": "biologic",
     "required_species": ["anthocyanins", "chlorophyll_a+b"],
     "bonus_species": ["carotenoids", "zeaxanthin_PRI"], "base_score": 0.75,
     "notes": "Anthocyanins + Chl: stress-protection pigmentation (UV/drought/cold/pathogen response)"},
    {"name": "phycoerythrin_algae",      "category": "biologic",
     "required_species": ["phycoerythrin", "phycoerythrin"],
     "bonus_species": ["chlorophyll_a"], "base_score": 0.88,
     "notes": "Phycoerythrin doublet 495+545nm: red algae or PE-type cyanobacteria in water"},
    {"name": "cyanobacteria_bloom",      "category": "biologic",
     "required_species": ["phycocyanin", "allophycocyanin", "chlorophyll_a+b"],
     "bonus_species": ["liquid_water_OH"], "base_score": 0.87,
     "notes": "Phycocyanin + allophycocyanin + Chl: cyanobacteria bloom (HAB indicator)"},
    {"name": "diatom_bloom",             "category": "biologic",
     "required_species": ["fucoxanthin", "chlorophyll_a"],
     "bonus_species": ["liquid_water_OH"], "base_score": 0.82,
     "notes": "Fucoxanthin + Chl-a: diatom-dominated phytoplankton bloom (spring bloom signal)"},

    # ── IRON MINERALOGY ──────────────────────────────────────────────────
    {"name": "ferrihydrite",             "category": "mineral",
     "required_species": ["ferrihydrite", "ferrihydrite"],
     "bonus_species": ["jarosite"], "base_score": 0.83,
     "notes": "Ferrihydrite 430+750nm pair; short-range-order Fe3+ oxyhydroxide; AMD and mine-waste indicator"},
    {"name": "lepidocrocite_soil",       "category": "mineral",
     "required_species": ["lepidocrocite", "lepidocrocite"],
     "base_score": 0.80,
     "notes": "Lepidocrocite γ-FeOOH 475+1000nm pair; waterlogged gleyed soils; redox-oscillation indicator"},
    {"name": "acid_mine_drainage",       "category": "geological",
     "required_species": ["ferrihydrite", "jarosite"],
     "bonus_species": ["szomolnokite_Fe2+", "Fe3+_goethite"], "base_score": 0.85,
     "notes": "AMD assemblage: ferrihydrite+jarosite±szomolnokite; low-pH Fe-oxidation drainage signature"},
    {"name": "siderite",                 "category": "mineral",
     "required_species": ["siderite_Fe2+", "siderite_CO3"],
     "base_score": 0.81,
     "notes": "Siderite FeCO3: broad Fe2+ ~1080nm + CO3 ~2340nm; anoxic/reducing sedimentary facies"},

    # ── REE EXPANSION ────────────────────────────────────────────────────
    {"name": "REE_erbium",               "category": "mineral",
     "required_species": ["Er3+", "Er3+"],
     "bonus_species": ["Sm3+_or_Er3+"], "base_score": 0.88,
     "notes": "Er3+ 4f multiplet ≥2 of {650, 975, 1530nm}; 1530nm telecom band is highly diagnostic"},
    {"name": "REE_dysprosium",           "category": "mineral",
     "required_species": ["Dy3+", "Dy3+"],
     "base_score": 0.87,
     "notes": "Dy3+ doublet 910+1260nm; HREE; carbonatite/granitic REE deposit indicator"},
    {"name": "REE_ytterbium",            "category": "mineral",
     "required_species": ["Yb3+"],
     "base_score": 0.72,
     "notes": "Yb3+ 975nm single narrow line; requires FWHM <15nm to distinguish from H2O atmospheric"},
    {"name": "REE_holmium",              "category": "mineral",
     "required_species": ["Ho3+", "Ho3+"],
     "base_score": 0.85,
     "notes": "Ho3+ doublet 450+537nm; paired narrow 4f features diagnostic of holmium"},

    # ── CLAY / PHYLLOSILICATE ────────────────────────────────────────────
    {"name": "pyrophyllite",             "category": "mineral",
     "required_species": ["Al-OH_pyrophyllite", "Al-OH_pyrophyllite"],
     "required_absent": ["Al-OH_clay"],
     "base_score": 0.86,
     "notes": "Pyrophyllite 2165+1390nm doublet; absent kaolinite 1410nm; epithermal acid-sulfate alteration mapping"},
    {"name": "gibbsite_laterite",        "category": "mineral",
     "required_species": ["gibbsite", "gibbsite"],
     "base_score": 0.85,
     "notes": "Gibbsite Al(OH)3 doublet 2263+2387nm; tropical laterite and deep weathering profile"},
    {"name": "boehmite_bauxite",         "category": "mineral",
     "required_species": ["boehmite", "boehmite"],
     "bonus_species": ["gibbsite"], "base_score": 0.83,
     "notes": "Boehmite AlOOH doublet 2175+2110nm; bauxite deposits; high-temperature laterite horizon"},
    {"name": "halloysite",               "category": "mineral",
     "required_species": ["halloysite", "liquid_water_OH"],
     "bonus_species": ["OH_water"], "base_score": 0.76,
     "notes": "Halloysite hydrated kaolinite: Al-OH 2205nm + interlayer water signature"},
    {"name": "talc",                     "category": "mineral",
     "required_species": ["talc_Mg-OH", "talc_Mg-OH"],
     "bonus_species": ["talc_OH"], "base_score": 0.86,
     "notes": "Talc Mg3Si4O10(OH)2: doublet 2315+2390nm diagnostic; ultramafic talc-carbonate alteration zones"},

    # ── SULPHATE MINERALS ────────────────────────────────────────────────
    {"name": "epsomite",                 "category": "mineral",
     "required_species": ["epsomite_H2O", "epsomite_MgSO4"],
     "base_score": 0.84,
     "notes": "Epsomite MgSO4·7H2O: water 1450+1940nm + Mg-SO4 2200nm; saline playa evaporite"},
    {"name": "kieserite",                "category": "mineral",
     "required_species": ["kieserite_H2O", "kieserite_MgSO4"],
     "base_score": 0.83,
     "notes": "Kieserite MgSO4·H2O: 1640nm monohydrate + 2220nm SO4; distinct from epsomite heptahydrate"},
    {"name": "szomolnokite",             "category": "mineral",
     "required_species": ["szomolnokite_Fe2+", "szomolnokite_H2O"],
     "base_score": 0.80,
     "notes": "Szomolnokite FeSO4·H2O: Fe2+ 1000nm + H2O 1490nm; AMD iron-sulfate mineral assemblage"},

    # ── GEOLOGICAL MINERALS ──────────────────────────────────────────────
    {"name": "prehnite",                 "category": "mineral",
     "required_species": ["prehnite_Al-OH", "prehnite_OH"],
     "base_score": 0.81,
     "notes": "Prehnite Ca2Al2Si3O10(OH)2: 2243+1490nm; low-grade metamorphic / prehnite-pumpellyite facies"},
    {"name": "schorl_tourmaline",        "category": "mineral",
     "required_species": ["schorl_B-OH", "Fe3+_charge_transfer"],
     "base_score": 0.78,
     "notes": "Schorl tourmaline: B-OH 2120nm + Fe3+ CT; metasediment, boron-rich greisen alteration"},
    {"name": "crocidolite_asbestos",     "category": "mineral",
     "required_species": ["crocidolite_Fe-OH", "crocidolite_Fe3+"],
     "base_score": 0.86,
     "notes": "HAZARD: crocidolite (blue asbestos) Fe-OH 2320nm + Fe3+ 440nm; Fe3+ VIS band absent in chrysotile"},
    {"name": "ankerite",                 "category": "mineral",
     "required_species": ["ankerite_Mg-OH", "ankerite_CO3"],
     "base_score": 0.80,
     "notes": "Ankerite Ca(Mg,Fe)(CO3)2: doublet 2320+2335nm; carbonate alteration halo in gold deposits"},

    # ── PLASTICS EXPANSION ───────────────────────────────────────────────
    {"name": "PET_plastic",              "category": "synthetic",
     "required_species": ["PET_CH", "PET_aromatic"],
     "base_score": 0.85,
     "notes": "PET polyester: 1730nm CH + 2170nm aromatic combination; bottles, textiles, litter mapping"},
    {"name": "PVC_plastic",              "category": "synthetic",
     "required_species": ["PVC_CH2", "PVC_CCl"],
     "base_score": 0.87,
     "notes": "PVC: 1730nm CH2 + 2030nm C-Cl combination; diagnostic C-Cl unique among common plastics"},
    {"name": "polystyrene",              "category": "synthetic",
     "required_species": ["polystyrene_CH_aromatic", "polystyrene_aromatic"],
     "base_score": 0.84,
     "notes": "Polystyrene: 1680nm aromatic CH + 2240nm ring combination; packaging, foam, litter mapping"},
    {"name": "rubber_polyisoprene",      "category": "synthetic",
     "required_species": ["rubber_C=C", "rubber_CH"],
     "base_score": 0.79,
     "notes": "Natural/synthetic rubber polyisoprene: C=C 1640nm + CH2 1720nm; tire debris mapping"},

    # ── URBAN MATERIALS ──────────────────────────────────────────────────
    {"name": "concrete_urban",           "category": "urban",
     "required_species": ["portlandite_CaOH", "carbonate_CO3"],
     "bonus_species": ["OH_water"], "base_score": 0.78,
     "notes": "Concrete/cement: portlandite Ca-OH 1460nm + carbonation CO3 ~2340nm; urban surface mapping"},
    {"name": "asphalt_bitumen",          "category": "urban",
     "required_species": ["CH2_lipid", "CH2_CH3_hydrocarbon", "coal_aromatic"],
     "required_absent": ["red_edge_inflection", "Al-OH"],
     "bonus_species": ["coal_CH2"], "base_score": 0.83,
     "notes": "Asphalt bitumen: CH2 1720nm + CH 2310nm + aromatic 670nm; absent vegetation/clay confirms road surface"},
]

# Build ordered class labels (0 = unknown, 1..N = composite rules)
CLASS_LABELS: list[str] = ["unknown"] + [r["name"] for r in COMPOSITE_RULES]
CLASS_CATEGORIES: list[str] = [""] + [r["category"] for r in COMPOSITE_RULES]


# ---------------------------------------------------------------------------
# Spectral interpretation engine (single-spectrum / point mode)
# ---------------------------------------------------------------------------

_CONTINUUM_MARGIN = 150.0  # nm half-window for continuum estimation


def _band_depth(feat_reflectance: float, feat_center: float,
                all_centers: list[float], all_refl: list[float]) -> float:
    """Band depth = 1 - R_band / R_continuum (local max in ±margin window)."""
    neighbors = [
        r for c, r in zip(all_centers, all_refl)
        if abs(c - feat_center) < _CONTINUUM_MARGIN and c != feat_center
    ]
    continuum = max(neighbors) if neighbors else 1.0
    if continuum <= 0:
        return 0.0
    return max(0.0, 1.0 - feat_reflectance / continuum)


def interpret_spectrum(
    triplets: list[tuple[float, float, float]],
    min_band_depth: float = 0.02,
) -> SpectralInterpretation:
    """Interpret a reflectance spectrum from absorption feature triplets.

    Parameters
    ----------
    triplets : list of (center_nm, fwhm_nm, reflectance)
    min_band_depth : minimum band depth to count as absorption

    Returns
    -------
    SpectralInterpretation with ranked material hypotheses.
    """
    flags: list[str] = []

    # Step 0: Build AbsorptionFeature objects
    raw: list[AbsorptionFeature] = []
    for idx, (c, fw, ref) in enumerate(triplets):
        if not (400 <= c <= 2500):
            flags.append(f"Band #{idx} center {c:.0f}nm out of 400–2500nm – skipped")
            continue
        if fw < 0.5:
            continue
        raw.append(AbsorptionFeature(
            center_nm=c, fwhm_nm=fw,
            reflectance=max(0.0, min(1.0, ref))
        ))
    raw.sort(key=lambda f: f.center_nm)

    # Step 1: Band depth from local continuum
    centers = [f.center_nm for f in raw]
    refls   = [f.reflectance for f in raw]
    for feat in raw:
        feat.depth = _band_depth(feat.reflectance, feat.center_nm, centers, refls)
        if feat.fwhm_nm < 10:
            feat.asymmetry = "sharp"
        elif feat.fwhm_nm < 50:
            feat.asymmetry = "moderate"
        else:
            feat.asymmetry = "broad"

    features = [f for f in raw if f.depth >= min_band_depth]
    if not features:
        flags.append("No features above min_band_depth. Returning empty result.")
        return SpectralInterpretation(
            features=raw, chromophores=[], hypotheses=[],
            dominant_class=None, atomic_summary=[], flag_notes=flags)

    # Step 2: Match features → DB species
    species_map: dict[str, list[tuple[AbsorptionFeature, dict]]] = {}
    for feat in features:
        for entry in SPECTRAL_FEATURES_DB:
            lo, hi = entry["range"]
            if lo <= feat.center_nm <= hi and feat.depth >= min_band_depth:
                species_map.setdefault(entry["species"], []).append((feat, entry))

    # Step 3: ChromophoreAssignment
    chromophores: list[ChromophoreAssignment] = []
    for sp, fe_list in species_map.items():
        matched = [fe[0] for fe in fe_list]
        n = len(matched)
        mean_depth = sum(f.depth for f in matched) / n
        depth_score = min(1.0, mean_depth / 0.3)
        multi_bonus = min(1.0, 0.5 + 0.25 * n)
        conf = round(depth_score * multi_bonus, 3)
        chromophores.append(ChromophoreAssignment(
            species=sp,
            mechanism=fe_list[0][1]["mechanism"],
            supporting_features=sorted(set(f.center_nm for f in matched)),
            confidence=conf,
            notes=fe_list[0][1]["notes"],
        ))
    chromophores.sort(key=lambda c: c.confidence, reverse=True)

    # Step 4: Composite rules → MaterialHypothesis
    detected = set(species_map.keys())
    hypotheses: list[MaterialHypothesis] = []

    for rule in COMPOSITE_RULES:
        req = rule["required_species"]
        req_counts: dict[str, int] = {}
        for sp in req:
            req_counts[sp] = req_counts.get(sp, 0) + 1
        if not all(
            sp in detected and len(species_map.get(sp, [])) >= cnt
            for sp, cnt in req_counts.items()
        ):
            continue
        if any(sp in detected for sp in rule.get("required_absent", [])):
            continue

        bonus_hit = sum(1 for sp in rule.get("bonus_species", []) if sp in detected)
        score = min(1.0, rule["base_score"] * (1.0 + 0.05 * bonus_hit))

        evidence = [
            f"{sp} @ {fe[0].center_nm:.0f}nm (depth={fe[0].depth:.3f})"
            for sp in set(req) for fe in species_map.get(sp, [])
        ]
        hypotheses.append(MaterialHypothesis(
            name=rule["name"], category=rule["category"],
            confidence=round(score, 3), evidence=evidence,
            atomic_species=list(set(req) | (set(rule.get("bonus_species", [])) & detected)),
        ))

    # Step 5: Atomic fallback for unmatched species
    matched_in_rules = {sp for h in hypotheses for sp in h.atomic_species}
    for ca in chromophores:
        if ca.species not in matched_in_rules and ca.confidence > 0.2:
            tags = next((e["tags"] for e in SPECTRAL_FEATURES_DB
                         if e["species"] == ca.species), [])
            hypotheses.append(MaterialHypothesis(
                name=f"atomic_{ca.species}",
                category=tags[0] if tags else "unknown",
                confidence=round(ca.confidence * 0.7, 3),
                evidence=[f"{ca.species} @ {w:.0f}nm" for w in ca.supporting_features],
                atomic_species=[ca.species],
            ))

    # Step 6: Atmospheric contamination flag
    atm = {"O2_A_band_atmospheric", "O2_B_band_atmospheric", "H2O_atmospheric"}
    if atm & detected:
        flags.append(
            f"Atmospheric features detected: {atm & detected}. "
            "Verify atmospheric correction before material interpretation."
        )

    # Step 7: SWIR overlap warning (2280–2390 nm)
    if any(2280 <= f.center_nm <= 2390 for f in features):
        ambig = [s for s in detected
                 if any(k in s for k in ("carbonate", "CH2_CH3", "Mg-OH", "lignin"))]
        if len(ambig) > 1:
            flags.append(
                f"Spectral overlap 2280–2390nm: {ambig}. "
                "Use FWHM (<30nm=organic; >80nm=mineral) to disambiguate."
            )

    hypotheses.sort(key=lambda h: h.confidence, reverse=True)
    dominant = None
    if hypotheses and hypotheses[0].confidence >= 0.60:
        dominant = hypotheses[0].name
    elif hypotheses:
        flags.append("No hypothesis reached 60% confidence.")

    return SpectralInterpretation(
        features=features, chromophores=chromophores, hypotheses=hypotheses,
        dominant_class=dominant,
        atomic_summary=sorted({sp for h in hypotheses for sp in h.atomic_species}),
        flag_notes=flags,
    )

# ---------------------------------------------------------------------------
# Vectorized per-pixel scoring (full-map mode)
# ---------------------------------------------------------------------------


def _compute_depths_numpy(cube: np.ndarray, wavelengths: np.ndarray,
                          margin: float = 150.0) -> np.ndarray:
    """Compute band depths for every pixel in the cube simultaneously.

    Parameters
    ----------
    cube        : shape (Z, pixels), float64, reflectance 0–1, NaN = nodata
    wavelengths : shape (Z,), band centre wavelengths in nm
    margin      : continuum half-window (nm)

    Returns
    -------
    depths : shape (Z, pixels), band depth = 1 - R / R_continuum, ≥ 0
    """
    Z = len(wavelengths)
    depths = np.zeros_like(cube)
    for i in range(Z):
        mask = np.abs(wavelengths - wavelengths[i]) < margin
        mask[i] = False
        if mask.any():
            local_max = np.nanmax(cube[mask], axis=0)
        else:
            local_max = np.ones(cube.shape[1])
        local_max = np.maximum(local_max, 1e-10)
        depths[i] = np.maximum(0.0, 1.0 - cube[i] / local_max)
    return depths


def _build_band_species_index(wavelengths: np.ndarray) -> dict[str, list[int]]:
    """For each DB species, list of band indices that fall within its range."""
    idx: dict[str, list[int]] = {}
    for entry in SPECTRAL_FEATURES_DB:
        sp = entry["species"]
        lo, hi = entry["range"]
        hits = [i for i, wl in enumerate(wavelengths) if lo <= wl <= hi]
        if hits:
            idx.setdefault(sp, []).extend(hits)
    # Deduplicate
    return {sp: sorted(set(v)) for sp, v in idx.items()}


def score_composites_numpy(depths: np.ndarray, wavelengths: np.ndarray,
                            min_depth: float = 0.02) -> np.ndarray:
    """Score all composite rules per pixel using vectorized numpy.

    Parameters
    ----------
    depths      : shape (Z, pixels), band depths from _compute_depths_numpy
    wavelengths : shape (Z,)
    min_depth   : minimum depth to count a species as detected

    Returns
    -------
    scores : shape (n_rules, pixels), score 0–1 per rule per pixel
    """
    n_pixels = depths.shape[1]
    sp_index = _build_band_species_index(wavelengths)
    n_rules = len(COMPOSITE_RULES)
    scores = np.zeros((n_rules, n_pixels), dtype=np.float32)

    # Pre-compute detection map: species → bool array (n_pixels,)
    def _species_detected(sp: str) -> np.ndarray:
        band_idxs = sp_index.get(sp, [])
        if not band_idxs:
            return np.zeros(n_pixels, dtype=bool)
        return np.any(depths[band_idxs] >= min_depth, axis=0)

    def _species_n_detections(sp: str) -> np.ndarray:
        band_idxs = sp_index.get(sp, [])
        if not band_idxs:
            return np.zeros(n_pixels, dtype=np.int32)
        return np.sum(depths[band_idxs] >= min_depth, axis=0).astype(np.int32)

    for r_idx, rule in enumerate(COMPOSITE_RULES):
        req = rule["required_species"]
        # Count how many times each species is required
        req_counts: dict[str, int] = {}
        for sp in req:
            req_counts[sp] = req_counts.get(sp, 0) + 1

        # All required species must be detected sufficient times
        all_met = np.ones(n_pixels, dtype=bool)
        for sp, cnt in req_counts.items():
            if cnt == 1:
                all_met &= _species_detected(sp)
            else:
                all_met &= (_species_n_detections(sp) >= cnt)

        # Required absent species
        for sp in rule.get("required_absent", []):
            all_met &= ~_species_detected(sp)

        if not np.any(all_met):
            continue

        # Bonus species
        bonus_count = np.zeros(n_pixels, dtype=np.float32)
        for sp in rule.get("bonus_species", []):
            bonus_count += _species_detected(sp).astype(np.float32)

        bonus_factor = 1.0 + 0.05 * bonus_count
        score = np.minimum(1.0, rule["base_score"] * bonus_factor).astype(np.float32)
        scores[r_idx] = np.where(all_met, score, 0.0)

    return scores

# ---------------------------------------------------------------------------
# Pretty printer (point mode)
# ---------------------------------------------------------------------------


def print_interpretation(result: SpectralInterpretation, top_n: int = 8,
                          verbose: bool = False) -> None:
    hr = "─" * 70
    gs.message(hr)
    gs.message("  SPECTRAL INTERPRETATION REPORT")
    gs.message(hr)
    gs.message(f"  Dominant class : {result.dominant_class or '(none)'}")
    gs.message(f"  Active features: {len(result.features)}")
    gs.message(f"  Chromophores   : {len(result.chromophores)}")
    gs.message(f"  Hypotheses     : {len(result.hypotheses)}")

    for fl in result.flag_notes:
        gs.warning(fl)

    gs.message(f"\n  RANKED HYPOTHESES (top {top_n})")
    for i, h in enumerate(result.hypotheses[:top_n]):
        bar = "█" * int(h.confidence * 20)
        gs.message(f"  [{i+1:2d}] {h.name:<35} {h.category:<12} "
                   f"{h.confidence:.2f}  {bar}")
        if verbose:
            for ev in h.evidence[:3]:
                gs.message(f"         ↳ {ev}")

    gs.message(f"\n  CHROMOPHORE DETECTIONS")
    for ca in result.chromophores[:12]:
        wls = ", ".join(f"{w:.0f}nm" for w in ca.supporting_features)
        gs.message(f"  {ca.species:<35} conf={ca.confidence:.2f}  [{wls}]")

    gs.message(hr)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(options, flags):
    input_map    = options['input']
    output_map   = options['output']
    conf_map     = options['confidence']
    min_depth    = float(options['min_band_depth'])
    coordinates  = options.get('coordinates', '')
    only_valid   = flags['n']
    info_mode    = flags['i']
    point_mode   = flags['p']
    verbose      = flags['v']

    # Collect band metadata
    bands = get_band_info(input_map, only_valid=only_valid)
    wavelengths = np.array([b['wavelength'] for b in bands], dtype=np.float64)
    fwhms       = np.array([b['fwhm']       for b in bands], dtype=np.float64)

    # ── Info mode ──────────────────────────────────────────────────────────
    if info_mode:
        gs.message(f"3D raster : {input_map}")
        gs.message(f"Bands     : {len(bands)}")
        gs.message(f"WL range  : {wavelengths[0]:.1f} – {wavelengths[-1]:.1f} nm")
        gs.message(f"{'Band':>5}  {'WL (nm)':>10}  {'FWHM':>8}  {'Valid':>5}")
        for b in bands:
            gs.message(f"  {b['band']:3d}  {b['wavelength']:10.2f}  "
                       f"{b['fwhm']:8.2f}  {'yes' if b['valid'] else 'no':>5}")
        return 0

    # ── Point mode ─────────────────────────────────────────────────────────
    if point_mode:
        if not coordinates:
            gs.fatal("Point mode requires coordinates=east,north")
        try:
            east, north = [float(v) for v in coordinates.split(',')]
        except ValueError:
            gs.fatal("coordinates must be 'east,north' (two floats)")

        gs.message(f"Extracting spectrum at E={east} N={north}")
        gs.verbose("Loading all bands via fast Z-slice extraction…")
        spectrum: list[tuple[float, float, float]] = []
        for b in bands:
            tmp = extract_band(input_map, b['band'])
            val_str = gs.read_command(
                'r.what', map=tmp, coordinates=f"{east},{north}", quiet=True
            ).split('|')[-1].strip()
            try:
                refl = float(val_str)
            except ValueError:
                refl = float('nan')
            spectrum.append((b['wavelength'], b['fwhm'], refl))

        # Filter out nodata
        spectrum = [(c, fw, r) for c, fw, r in spectrum if not (r != r)]  # NaN check
        if not spectrum:
            gs.fatal("All band values are nodata at the specified coordinates.")

        result = interpret_spectrum(spectrum, min_band_depth=min_depth)
        print_interpretation(result, verbose=verbose)
        return 0

    # ── Full-map mode ───────────────────────────────────────────────────────
    if not output_map and not conf_map:
        gs.fatal("Specify output= and/or confidence= for full-map mode.")

    info = gs.raster3d_info(input_map)
    nrows = int(info['rows'])
    ncols = int(info['cols'])
    n_pixels = nrows * ncols
    Z = len(bands)

    gs.message(f"Loading {Z} bands ({nrows}×{ncols} pixels) via fast Z-slice…")
    cube = np.full((Z, n_pixels), np.nan, dtype=np.float64)

    for z_idx, b in enumerate(bands):
        gs.percent(z_idx, Z, 2)
        tmp = extract_band(input_map, b['band'])
        # Read 2D raster into numpy array via r.out.bin or grass.pygrass
        try:
            from grass.pygrass.raster import RasterRow
            with RasterRow(tmp) as rmap:
                arr = np.array(rmap, dtype=np.float64).ravel()
        except Exception:
            # Fallback: write ASCII and read back
            ascii_out = gs.tempfile()
            gs.run_command('r.out.ascii', input=tmp, output=ascii_out,
                           null_value='nan', quiet=True)
            arr = np.loadtxt(ascii_out, dtype=np.float64).ravel()
            os.unlink(ascii_out)
        if arr.size == n_pixels:
            cube[z_idx] = arr
        else:
            gs.warning(f"Band {b['band']}: unexpected array size {arr.size}, expected {n_pixels}")

    gs.percent(Z, Z, 2)
    gs.message("Computing band depths…")
    depths = _compute_depths_numpy(cube, wavelengths, margin=_CONTINUUM_MARGIN)

    gs.message("Scoring composite rules per pixel…")
    scores = score_composites_numpy(depths, wavelengths, min_depth=min_depth)
    # scores: (n_rules, n_pixels)

    # Class map: argmax+1 (1-based), 0 where all scores == 0
    best_idx = np.argmax(scores, axis=0)          # (n_pixels,)  0-based rule idx
    best_score = scores[best_idx, np.arange(n_pixels)]

    class_map = (best_idx + 1).reshape(nrows, ncols).astype(np.int32)
    class_map[best_score.reshape(nrows, ncols) < 0.1] = 0

    conf_arr = best_score.reshape(nrows, ncols).astype(np.float32)
    # Set nodata pixels to 0
    nodata_mask = np.all(np.isnan(cube), axis=0).reshape(nrows, ncols)
    class_map[nodata_mask] = 0
    conf_arr[nodata_mask] = np.nan

    # Write output rasters
    def _write_raster(name, array, rtype, null_val):
        ascii_tmp = gs.tempfile()
        if rtype == 'int':
            np.savetxt(ascii_tmp, array, fmt='%d', delimiter=' ')
        else:
            np.savetxt(ascii_tmp, array, fmt='%.4f', delimiter=' ',
                       header='', comments='')
        gs.run_command('r.in.ascii', input=ascii_tmp, output=name,
                       null_value=str(null_val), overwrite=True, quiet=True)
        os.unlink(ascii_tmp)

    if output_map:
        gs.message(f"Writing class map → {output_map}")
        _write_raster(output_map, class_map, 'int', 0)

        # Write category labels
        cats = "\n".join(f"{i}:{CLASS_LABELS[i]}" for i in range(len(CLASS_LABELS)))
        gs.write_command('r.category', map=output_map,
                         rules='-', separator=':', stdin=cats, quiet=True)

        gs.run_command('r.support', map=output_map,
                       title="Spectral interpretation – dominant class",
                       units="class code",
                       description=(
                           f"i.hyper.spectroscopy class map. "
                           f"0=unknown, 1..{len(COMPOSITE_RULES)}=composite rules. "
                           f"Input: {input_map}"))

    if conf_map:
        gs.message(f"Writing confidence map → {conf_map}")
        _write_raster(conf_map, conf_arr, 'float', 'nan')
        gs.run_command('r.support', map=conf_map,
                       title="Spectral interpretation – confidence",
                       units="score 0-1",
                       description=(
                           f"i.hyper.spectroscopy confidence. "
                           f"Input: {input_map}"))
        gs.run_command('r.colors', map=conf_map, color='plasma', quiet=True)

    # Summary statistics
    if output_map:
        try:
            unique = np.unique(class_map[~nodata_mask])
            gs.message("Class distribution:")
            for code in unique:
                if code == 0:
                    continue
                cnt = np.sum(class_map == code)
                label = CLASS_LABELS[code] if code < len(CLASS_LABELS) else "?"
                gs.message(f"  {code:3d}  {label:<35}  {cnt:8d} px")
        except Exception:
            pass

    return 0


if __name__ == "__main__":
    options, flags = gs.parser()
    sys.exit(main(options, flags))
