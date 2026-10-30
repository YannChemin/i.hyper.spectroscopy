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
# % key: validity
# % type: string
# % required: no
# % label: Output per-pixel validity raster name
# % description: 1=classified, 2=unclassified, 3=low signal (below min_reflectance), 4=saturated (above max_reflectance); NULL where the cube has no data. Separates "no class could be supported here" from "a class was ruled out"
# % guisection: Output
# %end

# %option
# % key: reflectance_scale
# % type: double
# % required: no
# % answer: 0
# % label: Divisor that converts stored values to 0-1 reflectance (0 = detect from the data)
# % description: Surface-reflectance products are often stored scaled by 1000 or 10000. Band depth is a ratio and does not care, but the absolute thresholds (min_reflectance, max_reflectance, sand_albedo_threshold) do. Detected from the cube's own value range when 0, and reported
# % guisection: Processing
# %end

# %option
# % key: min_reflectance
# % type: double
# % required: no
# % answer: 0.02
# % label: Reflectance below which a band depth is treated as unmeasurable
# % description: A species whose bands are all darker than this is reported as indeterminate rather than detected or absent, and no rule that needs it can fire. Band depth is a ratio, so on a dark surface a large depth can correspond to a reflectance difference smaller than the sensor noise. Set to 0 to disable
# % guisection: Processing
# %end

# %option
# % key: max_reflectance
# % type: double
# % required: no
# % answer: 1.2
# % label: Reflectance above which a pixel is treated as saturated
# % description: Pixels with any band above this are flagged saturated in the validity map and left unclassified, instead of being reported as an ordinary no-match. Set very high (e.g. 100) to disable
# % guisection: Processing
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

# %option
# % key: ndre_threshold
# % type: double
# % required: no
# % answer: 0.15
# % description: Minimum NDRE (NIR vs red-edge) required to allow chlorophyll-dependent classes (vegetation, algae); screens out false positives on bare soil/rock. Set very negative (e.g. -1) to disable.
# % guisection: Processing
# %end

# %option
# % key: sand_albedo_threshold
# % type: double
# % required: no
# % answer: 0.25
# % description: Minimum broadband 500-900nm reflectance for the quartz/feldspar sand fallback (quartz and feldspar have no diagnostic VNIR-SWIR absorption feature, so they are identified by elimination + brightness rather than a positive match). Set to 1.0 or higher to disable.
# % guisection: Processing
# %end

# %option
# % key: species_percentile
# % type: double
# % required: no
# % answer: 90
# % description: Full-map mode only. A species counts as detected only where its band depth is ALSO above this percentile of its own depth distribution across the scene (in addition to min_band_depth), so a broad slope common to the whole scene (e.g. iron-oxide/clay desert soil) does not register as a "feature" everywhere. Set to 0 to disable (absolute thresholds only).
# % guisection: Processing
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

# ── ras3d standalone detection ────────────────────────────────────────────────
import os as _os
_RAS3D = False
if not _os.environ.get('GISBASE'):
    try:
        import importlib.util as _ilu
        if _ilu.find_spec('ras3d') and _ilu.find_spec('ras3d_grass_shim'):
            from ras3d_grass_shim import install as _r3_install
            _r3_install()
            _RAS3D = True
    except Exception:
        pass
# ─────────────────────────────────────────────────────────────────────────────

import sys
import os
import atexit
import warnings
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
# 3D raster Z-slice extraction
# ---------------------------------------------------------------------------


def extract_band(raster3d: str, band_num: int) -> str:
    """Extract band_num (1-based) from raster3d into a temporary 2D raster.

    Uses the standard r3.to.rast module, windowing the 3D region to a single
    depth slice. (Previously used a custom Rast3d_extract_z_slice() ctypes
    fast-path, but that corrupts memory -- malloc(): mismatching
    next->prev_size -- when the active region is a sub-window of the
    raster's native footprint rather than its full extent; r3.to.rast has no
    such issue and self-reconciles any 2D/3D region mismatch.)

    Returns the temporary 2D raster name (registered for cleanup at exit).
    """
    if _RAS3D:
        import ras3d as _r3
        import ras3d_write as _r3w
        from ras3d_grass_shim import get_band_cache as _gbc
        _h = _r3.open_cube(raster3d)
        _arr = _r3.get_band(_h, band_num - 1)
        _r3.close_cube(_h)
        tmp_name = f"tmp_spectroscopy_{raster3d.replace('@','_').replace('#','_').replace('.','_')}_{band_num}"
        _r3w.write_raster2d(_r3w.outpath(tmp_name), _arr, None)
        _gbc()[tmp_name] = _arr
        _TMP_RASTERS.append(tmp_name)
        return tmp_name

    z = band_num - 1  # 1-based → 0-based
    base = raster3d.replace('@', '_').replace('#', '_').replace('.', '_')
    tmp_prefix = f"tmp_spectroscopy_{base}_{band_num}"
    tmp_name = f"{tmp_prefix}_00001"  # r3.to.rast's single-slice output suffix

    saved = gs.parse_command('g.region', flags='3g')
    gs.run_command('g.region', b=z, t=z + 1, tbres=1)
    try:
        gs.run_command('r3.to.rast', input=raster3d, output=tmp_prefix,
                       overwrite=True, quiet=True)
    finally:
        gs.run_command('g.region', b=saved['b'], t=saved['t'],
                       tbres=saved['tbres'])

    _TMP_RASTERS.append(tmp_name)
    return tmp_name

# ---------------------------------------------------------------------------
# Band metadata
# ---------------------------------------------------------------------------


def _spectrum_via_r3what(raster3d: str, east: float,
                         north: float) -> Optional[list[float]]:
    """Read one pixel's whole depth column with a single r3.what call.

    Returns one value per depth of the cube, indexed by depth, so the caller
    reads band *n* as ``values[n - 1]``; None if r3.what is unavailable or
    returns an unexpected number of values, so the caller can fall back to
    slicing.

    r3.what reads only the voxels asked for, so a full spectrum comes back in
    well under a second. Slicing the cube into one 2D raster per band to query
    a single pixel costs minutes on a large cube, because every band pays for
    the whole scene.

    r3.what reports one value per depth of the computational region, not of
    the map, so the region's depth window is aligned to the cube's own stack
    first. That happens inside a temporary region, which leaves the caller's
    region untouched.
    """
    from grass.exceptions import CalledModuleError

    try:
        depths = int(gs.raster3d_info(raster3d)["depths"])
    except (CalledModuleError, KeyError, ValueError):
        return None
    try:
        gs.use_temp_region()
    except (CalledModuleError, OSError):
        return None
    try:
        gs.run_command("g.region", raster_3d=raster3d, quiet=True)
        out = gs.read_command("r3.what", input=raster3d,
                              coordinates=f"{east},{north}",
                              null_value="nan", quiet=True)
    except CalledModuleError:
        return None
    finally:
        gs.del_temp_region()
    line = next((ln for ln in out.splitlines() if ln.strip()), "")
    # east|north|depth 1|depth 2|... in the cube's own depth order.
    parts = line.split("|")[2:]
    if len(parts) != depths:
        gs.verbose(
            f"r3.what returned {len(parts)} values for a {depths}-depth cube; "
            "falling back to 2D slice extraction"
        )
        return None
    values = []
    for v in parts:
        try:
            values.append(float(v))
        except ValueError:
            values.append(float("nan"))
    return values


def _load_hyper_json_bands(raster3d: str) -> list[dict]:
    """Read wavelength/fwhm/validity from i.hyper.import's JSON sidecar.

    i.hyper.import (HyperMetadata) stores band metadata at
    $MAPSET/grid3/<mapname>/hyper.json rather than in r3.support history,
    so that must be checked before falling back to r3.info/r.support.
    """
    import json as _json

    name, mapset = (raster3d.split('@', 1) if '@' in raster3d
                     else (raster3d, None))
    try:
        env = gs.gisenv()
        mapset = mapset or env['MAPSET']
        path = os.path.join(env['GISDBASE'], env['LOCATION_NAME'], mapset,
                            'grid3', name, 'hyper.json')
    except Exception:
        return []

    if not os.path.isfile(path):
        return []

    with open(path) as _fj:
        data = _json.load(_fj)

    b = data.get('bands') or {}
    wavelengths = b.get('wavelength')
    if not wavelengths:
        return []
    fwhms = b.get('fwhm') or []
    valids = b.get('validity') or []

    bands = []
    for i, wl in enumerate(wavelengths):
        bands.append({
            'band': i + 1,
            'wavelength': float(wl),
            'fwhm': float(fwhms[i]) if i < len(fwhms) else 10.0,
            'valid': bool(valids[i]) if i < len(valids) else True,
        })
    return bands


def get_band_info(raster3d: str, only_valid: bool = False) -> list[dict]:
    """Return sorted list of {band, wavelength, fwhm, valid} dicts.

    Primary source is i.hyper.import's grid3/<map>/hyper.json sidecar.
    Falls back to r3.info -h history (format: 'Band N: W nm, FWHM: F nm'),
    then to r.support metadata per band if history is missing.
    """
    if _RAS3D:
        import json as _json
        _sidecar = raster3d + '.wl.json'
        if os.path.isfile(_sidecar):
            with open(_sidecar) as _fj:
                _wl_data = _json.load(_fj)
            _bands = []
            for _entry in _wl_data:
                _bands.append({
                    'band': int(_entry['band']),
                    'wavelength': float(_entry['wavelength']),
                    'fwhm': float(_entry.get('fwhm', 10.0)),
                    'valid': bool(_entry.get('valid', True)),
                })
            _bands.sort(key=lambda b: b['wavelength'])
            if only_valid:
                _bands = [b for b in _bands if b['valid']]
            if _bands:
                return _bands

    json_bands = _load_hyper_json_bands(raster3d)
    if json_bands:
        json_bands.sort(key=lambda b: b['wavelength'])
        if only_valid:
            json_bands = [b for b in json_bands if b['valid']]
            if not json_bands:
                gs.fatal("No valid bands found (all marked valid=0).")
        return json_bands

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

    # Camouflage material. The 1500-1800 nm entries come from USGS
    # splib07 laboratory fabrics; the 2122/2245 nm pair and the water
    # exclusion were measured on the C2RD targets themselves.
    {"center": 1202, "range": (1192, 1222), "species": "coated_fabric_CH_1202",
     "mechanism": "overtone_CH",
     "tags": ["synthetic", "polymer", "fabric", "camouflage"],
     "notes": "C-H second overtone of the fabric coating; present in every synthetic sample"},
    {"center": 1414, "range": (1402, 1428), "species": "coated_fabric_NH_OH_1414",
     "mechanism": "overtone_NH_OH",
     "tags": ["synthetic", "polymer", "fabric", "camouflage"],
     "notes": "N-H / O-H first overtone; in all 8 synthetic samples, 1411-1419 nm"},
    {"center": 1535, "range": (1526, 1545), "species": "polyamide_NH_1535",
     "mechanism": "overtone_NH",
     "tags": ["synthetic", "polymer", "fabric", "camouflage"],
     "notes": "Polyamide N-H first overtone; nylon-specific, absent from cellulose"},
    {"center": 1573, "range": (1562, 1584), "species": "polyamide_NH_1573",
     "mechanism": "overtone_NH",
     "tags": ["synthetic", "polymer", "fabric", "camouflage"],
     "notes": "Second polyamide N-H component, always paired with the 1535 nm band"},
    {"center": 1716, "range": (1700, 1740), "species": "polymer_CH_overtone_1716",
     "mechanism": "overtone_CH",
     "tags": ["synthetic", "polymer", "fabric", "camouflage"],
     "notes": "C-H first overtone, the strongest synthetic-fabric band (depth 0.16-0.49)"},
    {"center": 1757, "range": (1742, 1772), "species": "polymer_CH_overtone_1755",
     "mechanism": "overtone_CH",
     "tags": ["synthetic", "polymer", "fabric", "camouflage"],
     "notes": "C-H first overtone partner; the 1716+1755 doublet is the branch's anchor"},
    {"center": 2049, "range": (2036, 2064), "species": "polyamide_amide_2049",
     "mechanism": "combination_amide",
     "tags": ["synthetic", "polymer", "fabric", "camouflage"],
     "notes": "Amide combination band of polyamide fabrics; depth up to 0.36"},
    {"center": 2183, "range": (2170, 2198), "species": "polymer_CH_2183",
     "mechanism": "combination_CH",
     "tags": ["synthetic", "polymer", "fabric", "camouflage"],
     "notes": "C-H combination of coated polyamide fabric"},
    {"center": 2294, "range": (2282, 2308), "species": "polymer_CH_2294",
     "mechanism": "combination_CH",
     "tags": ["synthetic", "polymer", "fabric", "camouflage"],
     "notes": "C-H combination; deepest SWIR band of the coated fabrics"},
    {"center": 2122, "range": (2105, 2140), "species": "polymer_sheet_2122",
     "mechanism": "combination_CH",
     "tags": ["synthetic", "polymer", "fabric", "camouflage"],
     "notes": "Field-measured band of the C2RD polymer sheeting"},
    {"center": 2245, "range": (2228, 2262), "species": "polymer_sheet_2245",
     "mechanism": "combination_CH",
     "tags": ["synthetic", "polymer", "fabric", "camouflage"],
     "notes": "Field-measured partner band; with 2122 nm it gives 11/12 recall"},
    # Open-water exclusion. Used only to veto synthetic-material rules, never
    # as a positive match, so every band carries its own min_depth: the OH
    # overtones differ in strength by a factor of two and a single species-
    # level threshold cannot serve all four. The 970 and 1200 nm floors are
    # measured on the C2RD Pourrieres scene (pool and pond water reach
    # 0.80-0.89 and 0.75-0.78; the deepest non-water surface in that scene,
    # a black polymer sheet, reaches 0.35 and 0.19). The 1450 nm band is
    # deliberately set out of reach: leaf water alone reaches 0.54 there, so
    # it cannot separate water from vegetation. The 1940 nm floor is carried
    # over from the laboratory branch and is not validated on field data,
    # because no sensor used here samples that window.
    {"center": 970, "range": (950, 990), "species": "water_body_OH_exclusion",
     "mechanism": "2nd_overtone_OH_stretch", "min_depth": 0.60,
     "tags": ["water", "exclusion"],
     "notes": "Open water 0.80-0.89 in the C2RD scene; deepest non-water surface 0.35"},
    {"center": 1200, "range": (1180, 1220), "species": "water_body_OH_exclusion",
     "mechanism": "combination_bend+stretch", "min_depth": 0.55,
     "tags": ["water", "exclusion"],
     "notes": "Open water 0.75-0.78 in the C2RD scene; deepest non-water surface 0.31"},
    {"center": 1450, "range": (1400, 1500), "species": "water_body_OH_exclusion",
     "mechanism": "1st_overtone_OH_stretch", "min_depth": 0.80,
     "tags": ["water", "exclusion"],
     "notes": "Inert by design: leaf water reaches 0.54 here, so 1450 nm cannot veto on its own"},
    {"center": 1940, "range": (1900, 1980), "species": "water_body_OH_exclusion",
     "mechanism": "combination_OH", "min_depth": 0.50,
     "tags": ["water", "exclusion"],
     "notes": "Laboratory value, not validated on field data: no sensor used here samples 1900-1980 nm"},

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

    # ── CAMOUFLAGE ────────────────────────────────────────────────────────
    {"name": "camouflage_fabric_cut_sheet", "category": "synthetic",
     "required_species": ["polymer_CH_overtone_1716", "polymer_CH_overtone_1755"],
     "min_depth_overrides": {"polymer_CH_overtone_1716": 0.12, "polymer_CH_overtone_1755": 0.1, "polyamide_NH_1535": 0.05, "polyamide_amide_2049": 0.05},
     "required_absent": ["water_body_OH_exclusion"],
     "bonus_species": ["coated_fabric_CH_1202", "coated_fabric_NH_OH_1414", "polymer_CH_2294"],
     "base_score": 0.78,
     "notes": "Die-cut sheet camouflage garnish: one coated fabric punched with leaf-shaped voids, not knotted netting"},
    {"name": "camouflage_fabric_polyamide", "category": "synthetic",
     "required_species": ["polymer_CH_overtone_1716", "polyamide_NH_1535", "polyamide_amide_2049"],
     "min_depth_overrides": {"polymer_CH_overtone_1716": 0.12, "polymer_CH_overtone_1755": 0.1, "polyamide_NH_1535": 0.05, "polyamide_amide_2049": 0.05},
     "required_absent": ["water_body_OH_exclusion"],
     "bonus_species": ["polyamide_NH_1573", "polymer_CH_2183", "polymer_CH_overtone_1755"],
     "base_score": 0.82,
     "notes": "Polyamide (nylon) camouflage fabric or webbing; amide bands separate it from polyolefin sheeting"},
    {"name": "camouflage_fabric_coated_polyolefin", "category": "synthetic",
     "required_species": ["polymer_CH_overtone_1716", "polymer_CH", "CH_aliphatic"],
     "min_depth_overrides": {"polymer_CH_overtone_1716": 0.12, "polymer_CH_overtone_1755": 0.1, "polyamide_NH_1535": 0.05, "polyamide_amide_2049": 0.05},
     "required_absent": ["water_body_OH_exclusion"],
     "bonus_species": ["coated_fabric_CH_1202", "coated_fabric_NH_OH_1414", "polymer_CH_overtone_1755"],
     "base_score": 0.8,
     "notes": "PVC- or polyolefin-coated woven camouflage sheet; very deep 2310/2350 nm bands, no amide"},
    {"name": "camouflage_garnish_over_vegetation", "category": "synthetic",
     "required_species": ["polymer_CH_overtone_1716", "chlorophyll_a+b"],
     "min_depth_overrides": {"polymer_CH_overtone_1716": 0.12, "polymer_CH_overtone_1755": 0.1, "polyamide_NH_1535": 0.05, "polyamide_amide_2049": 0.05},
     "required_absent": ["water_body_OH_exclusion"],
     "bonus_species": ["polyamide_NH_1535", "red_edge_inflection", "polymer_CH_overtone_1755"],
     "base_score": 0.62,
     "notes": "Cut-sheet garnish laid over vegetation: fabric C-H bands mixed with chlorophyll seen through the voids. The intended operational signature"},
    {"name": "camouflage_polymer_sheet_field", "category": "synthetic",
     "required_species": ["polymer_sheet_2122", "polymer_sheet_2245"],
     "min_depth_overrides": {"polymer_sheet_2122": 0.12, "polymer_sheet_2245": 0.1},
     "required_absent": ["water_body_OH_exclusion"],
     "bonus_species": ["polymer_CH_overtone_1716", "polymer_CH_2183"],
     "base_score": 0.75,
     "notes": "Polymer sheeting, derived from the C2RD Pourrieres targets themselves. Detects 11 of 12 sheets; the other hits are a bin liner, a transparent bag and a windsock, which are also polymer. Open water is excluded by its own deep 1940 nm band"},
    {"name": "natural_fibre_fabric", "category": "organic",
     "required_species": ["cellulose_OH", "cellulose"],
     "required_absent": ["polymer_CH_overtone_1716"],
     "bonus_species": ["lignin"],
     "base_score": 0.7,
     "notes": "Cotton, hessian or canvas. Separated from synthetic camouflage by the absence of the 1716 nm C-H overtone"},

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
     "required_species": ["ferrihydrite", "ferrihydrite", "jarosite"],
     "min_depth_overrides": {"jarosite": 0.08},
     "bonus_species": ["szomolnokite_Fe2+", "Fe3+_goethite"], "base_score": 0.90,
     "notes": "AMD assemblage: ferrihydrite (both 430+750nm bands)+jarosite±szomolnokite; "
              "low-pH Fe-oxidation drainage signature. Requiring only one ferrihydrite band "
              "and the default 0.02 jarosite depth let ordinary hematite/goethite (common, "
              "unrelated) false-positive via their shared generic blue-edge curvature and "
              "incidental shallow ~2265nm dips; tightened against USGS splib07 reference "
              "spectra (real jarosite/ferrihydrite: 0.11-0.17 depth @2265nm vs hematite/"
              "goethite's incidental 0.03-0.06). base_score raised to 0.90 (from 0.85, "
              "above the 'ferrihydrite' rule's max bonus-boosted 0.8715) so that when "
              "AMD's own stricter jarosite requirement is met, this more specific "
              "diagnosis deterministically outranks the generic ferrihydrite-alone label "
              "for the same pixel, instead of losing to it on an incidental bonus factor."},
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
# NDRE pre-screen for chlorophyll-dependent classes
# ---------------------------------------------------------------------------
#
# The "chlorophyll_a+b" species is matched by *any* shallow (>=min_band_depth)
# local reflectance dip in the 415-445nm or 670-695nm windows, and
# "red_edge_inflection" likewise by any dip in 700-730nm. Bare mineral/soil
# surfaces routinely show shallow local concavities in those same windows
# from unrelated causes (Fe3+ charge-transfer edge curvature, clay/carbonate
# combination-band shoulders, sensor/atmospheric-correction ripple), with
# nothing in the rule itself checking the actual NIR-vs-red-edge reflectance
# contrast that is what actually distinguishes live vegetation from bare
# ground. Combined with green_vegetation's base_score=0.90 (the highest of
# any rule), this lets bare desert pixels win the argmax purely on generic
# shallow-depth coincidences. NDRE (Normalized Difference Red Edge) is used
# here as an independent, physically-grounded pre-screen: pixels below the
# threshold are barred from winning any rule that requires "chlorophyll_a+b"
# (vegetation and chlorophyll-bearing algae/bloom classes), regardless of
# what the band-depth matching found.

_NDRE_RED_EDGE_NM = 710.0   # matches red_edge_inflection center (700-730nm)
_NDRE_NIR_NM = 865.0        # NIR shoulder, clear of major water-vapor bands

_CHLOROPHYLL_GATED_RULE_NAMES: set[str] = {
    r["name"] for r in COMPOSITE_RULES if "chlorophyll_a+b" in r["required_species"]
}
_CHLOROPHYLL_GATED_RULE_INDICES: list[int] = [
    i for i, r in enumerate(COMPOSITE_RULES) if r["name"] in _CHLOROPHYLL_GATED_RULE_NAMES
]


def _nearest_band_index(wavelengths: np.ndarray, target_nm: float) -> int:
    return int(np.argmin(np.abs(wavelengths - target_nm)))


def _safe_ndre(nir: np.ndarray, red_edge: np.ndarray) -> np.ndarray:
    denom = nir + red_edge
    with np.errstate(invalid='ignore', divide='ignore'):
        ndre = np.where(np.abs(denom) > 1e-6, (nir - red_edge) / denom, np.nan)
    return ndre

# ---------------------------------------------------------------------------
# Quartz/feldspar sand: elimination-based fallback (no diagnostic feature)
# ---------------------------------------------------------------------------
#
# Unlike every other class here, quartz and feldspar have no diagnostic
# electronic/vibrational absorption feature in 400-2500nm reflectance
# (confirmed against USGS splib07 reference spectra for quartz sand, albite,
# anorthite and microcline: in every case the only "feature" the same
# band-depth algorithm finds is the generic UV/blue absorption edge shared by
# nearly all minerals -- not something diagnostic). They are normally mapped
# via TIR emissivity (Christiansen feature, reststrahlen bands ~8-12um),
# outside this instrument's range. So instead of a positive species match,
# this class is assigned by elimination: no other rule scored above
# threshold, NDRE has already ruled out vegetation, and the pixel is bright
# (high broadband VNIR reflectance, as bare quartz/feldspar sand typically
# is). It is intentionally excluded from COMPOSITE_RULES / the normal
# per-pixel scoring loop and applied as a post-argmax fallback instead.

_SAND_BRIGHTNESS_RANGE_NM = (500.0, 900.0)
QUARTZOFELDSPATHIC_SAND_LABEL = "quartzofeldspathic_sand"
QUARTZOFELDSPATHIC_SAND_CONFIDENCE = 0.5  # modest: elimination, not a positive match


def _broadband_reflectance(cube_or_chunk: np.ndarray, bright_mask: np.ndarray) -> np.ndarray:
    # All-NaN columns (true no-data pixels) are expected here and produce NaN
    # output, which is fine (those pixels are masked out separately via
    # nodata_mask); suppress the resulting benign "Mean of empty slice"
    # RuntimeWarning, same rationale as _compute_depths_numpy's nanmax guard.
    with np.errstate(invalid='ignore'), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(cube_or_chunk[bright_mask], axis=0)


# Reserve a class number for the fallback, appended after (not part of)
# COMPOSITE_RULES so it never enters the normal per-pixel scoring loop.
CLASS_LABELS.append(QUARTZOFELDSPATHIC_SAND_LABEL)
CLASS_CATEGORIES.append("geological")
QUARTZOFELDSPATHIC_SAND_CLASS_IDX = len(CLASS_LABELS) - 1


# ---------------------------------------------------------------------------
# Spectral interpretation engine (single-spectrum / point mode)
# ---------------------------------------------------------------------------

_CONTINUUM_MARGIN = 150.0  # nm half-window for continuum estimation

# A continuum must be built from bands the sensor actually measured next to
# the feature. Most airborne and drone products drop the two water-vapour
# windows, leaving spectral gaps of 150-200 nm, and a plain "within +-150 nm"
# window then reaches across a gap and takes the continuum from the far side.
# That is not a shoulder of the feature at all, and it is how a water band
# depth of 0.64-0.70 came to be reported for a cube with no bands within
# 150 nm of 1940 nm. A gap is a spacing far larger than the sensor's own
# typical band spacing, so the test is relative: this keeps coarse
# multispectral sensors (where 60 nm between bands is normal) in one segment
# while still splitting a 5 nm imaging spectrometer at a 150 nm hole.
_GAP_MIN_NM = 50.0
_GAP_SPACING_FACTOR = 5.0

# Confidence calibration. A rule's base score says how diagnostic the rule
# is when it fires; on its own it says nothing about how well this particular
# pixel fits, which is why a strict run on the C2RD scene could report a mean
# confidence of 0.99 over what was a degenerate two-way split. The reported
# confidence is therefore the base score scaled by the depth margin of the
# rule's weakest required species: _EVIDENCE_DEPTH_SCALE is how much extra
# depth above the threshold counts as full evidence, and _EVIDENCE_FLOOR is
# the fraction of the base score a match earns when it only just clears every
# threshold. The scale is set from the module's own depth guide, where 0.10
# to 0.30 is a moderate feature and most composite rules activate there.
_EVIDENCE_DEPTH_SCALE = 0.20
_EVIDENCE_FLOOR = 0.55


def _spectral_segments(wavelengths: np.ndarray) -> np.ndarray:
    """Contiguous-coverage segment index for each band.

    Returns an array of the same length as *wavelengths*, numbering the
    contiguous blocks of coverage from 0 upwards; the index increases at each
    spectral gap. Bands in different segments are never used as each other's
    continuum. Assumes *wavelengths* is sorted ascending, as
    :func:`get_band_info` returns it; on unsorted input no gap is detected and
    the behaviour falls back to a single segment.
    """
    wl = np.asarray(wavelengths, dtype=float)
    if wl.size < 2:
        return np.zeros(wl.size, dtype=np.int32)
    spacing = np.diff(wl)
    threshold = max(_GAP_MIN_NM, _GAP_SPACING_FACTOR * float(np.median(spacing)))
    return np.concatenate(([0], np.cumsum(spacing > threshold))).astype(np.int32)


def _band_depth(feat_reflectance: float, feat_center: float,
                all_centers: list[float], all_refl: list[float]) -> float:
    """Band depth = 1 - R_band / R_continuum (local max in ±margin window).

    Only bands in the feature's own contiguous coverage segment can serve as
    its continuum, so a shoulder is never taken from across a spectral gap.
    A feature with no usable neighbour has no continuum and therefore no
    measurable depth: this returns 0.0 rather than assuming a continuum of
    1.0, which would have reported a depth of ``1 - R`` -- a large spurious
    value on any dark surface -- for an isolated band.
    """
    centers = np.asarray(all_centers, dtype=float)
    segments = _spectral_segments(centers)
    own = segments[int(np.argmin(np.abs(centers - feat_center)))]
    neighbors = [
        r for c, r, seg in zip(all_centers, all_refl, segments)
        if abs(c - feat_center) < _CONTINUUM_MARGIN and c != feat_center
        and seg == own
    ]
    if not neighbors:
        return 0.0
    continuum = max(neighbors)
    if continuum <= 0:
        return 0.0
    return max(0.0, 1.0 - feat_reflectance / continuum)


def interpret_spectrum(
    triplets: list[tuple[float, float, float]],
    min_band_depth: float = 0.02,
    ndre: Optional[float] = None,
    ndre_threshold: float = 0.15,
    wavelengths: Optional[np.ndarray] = None,
    min_reflectance: float = 0.0,
) -> SpectralInterpretation:
    """Interpret a reflectance spectrum from absorption feature triplets.

    Parameters
    ----------
    triplets : list of (center_nm, fwhm_nm, reflectance)
    min_band_depth : minimum band depth to count as absorption
    ndre : precomputed NDRE (NIR vs red-edge) for this pixel/spectrum, or
        None to skip the pre-screen entirely
    ndre_threshold : below this, chlorophyll-dependent composite rules
        (vegetation, algae blooms) are barred regardless of band-depth match

    Returns
    -------
    SpectralInterpretation with ranked material hypotheses.
    """
    flags: list[str] = []
    ndre_blocks_chlorophyll = (
        ndre is not None and not np.isnan(ndre) and ndre < ndre_threshold
    )
    if ndre_blocks_chlorophyll:
        flags.append(
            f"NDRE={ndre:.3f} < {ndre_threshold}: chlorophyll-dependent classes "
            f"({sorted(_CHLOROPHYLL_GATED_RULE_NAMES)}) screened out"
        )

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

    # Bands this spectrum actually carries per species, for the band-support
    # half of the confidence. Counted from the input, not from the matches,
    # so a species detected in one of eight sampled bands is not treated as
    # strongly as one detected in all eight.
    sp_band_counts: dict[str, list[float]] = {}
    for entry in SPECTRAL_FEATURES_DB:
        lo, hi = entry["range"]
        here = [c for c, _fw, _r in triplets if lo <= c <= hi]
        if here:
            sp_band_counts.setdefault(entry["species"], []).extend(here)

    # Rules the sensor cannot decide are dropped here too, for the same
    # reason as in the full-map path: an untestable exclusion would pass.
    unsupported = unsupported_rules(wavelengths) if wavelengths is not None else {}
    for name, reason in sorted(unsupported.items()):
        flags.append(f"Rule {name} disabled: {reason}")

    # Species too dark here for a band depth to carry information. Same test
    # as the full-map path: neither detected nor absent, so any rule needing
    # one is withheld rather than decided on noise.
    indeterminate: set[str] = set()
    if min_reflectance > 0:
        refl_by_species: dict[str, float] = {}
        for entry in SPECTRAL_FEATURES_DB:
            lo, hi = entry["range"]
            here = [r for c, _fw, r in triplets if lo <= c <= hi]
            if here:
                sp = entry["species"]
                refl_by_species[sp] = max(refl_by_species.get(sp, -np.inf), max(here))
        indeterminate = {sp for sp, r in refl_by_species.items()
                         if r < min_reflectance}
        for sp in sorted(indeterminate):
            flags.append(
                f"{sp} is indeterminate: its bands are below "
                f"{min_reflectance:g} reflectance, too dark for a band depth "
                "to be meaningful"
            )

    for rule in COMPOSITE_RULES:
        if rule["name"] in unsupported:
            continue
        if ndre_blocks_chlorophyll and rule["name"] in _CHLOROPHYLL_GATED_RULE_NAMES:
            continue
        req = rule["required_species"]
        overrides = rule.get("min_depth_overrides", {})
        req_counts: dict[str, int] = {}
        for sp in req:
            req_counts[sp] = req_counts.get(sp, 0) + 1
        def _threshold(sp: str, entry: dict, overrides: dict = overrides) -> float:
            """Deeper of the rule's override and the DB entry's own floor."""
            return max(overrides.get(sp, min_band_depth),
                       float(entry.get("min_depth", 0.0)))

        if not all(
            sum(1 for fe in species_map.get(sp, [])
                if fe[0].depth >= _threshold(sp, fe[1])) >= cnt
            for sp, cnt in req_counts.items()
        ):
            continue
        # Match the full-map path: an absent species honours its own depth
        # override, so a rule can require "no deep water band" rather than
        # "no water band at all".
        if any(
            any(fe[0].depth >= _threshold(sp, fe[1])
                for fe in species_map.get(sp, []))
            for sp in rule.get("required_absent", [])
        ):
            continue

        if indeterminate.intersection(req_counts) or indeterminate.intersection(
                rule.get("required_absent", [])):
            continue

        bonus_hit = sum(1 for sp in rule.get("bonus_species", []) if sp in detected)
        # Same calibration as the full-map path: the base score is the
        # ceiling and the weakest required species' depth margin decides how
        # much of it this spectrum earns.
        margins = []
        for sp in req_counts:
            matches = species_map.get(sp, [])
            best = max((fe[0].depth - _threshold(sp, fe[1])
                        for fe in matches), default=0.0)
            margin = min(1.0, max(0.0, best / _EVIDENCE_DEPTH_SCALE))
            # Band support: of the species' bands this sensor samples, how
            # many actually cleared the threshold here.
            available = len(sp_band_counts.get(sp, ())) or 1
            cleared = sum(1 for fe in matches
                          if fe[0].depth >= _threshold(sp, fe[1]))
            support = min(1.0, cleared / available)
            margins.append((margin * support) ** 0.5)
        evidence = (
            float(np.exp(np.mean(np.log(np.maximum(margins, 1e-6)))))
            if margins else 1.0
        )
        strength = _EVIDENCE_FLOOR + (1.0 - _EVIDENCE_FLOOR) * evidence
        score = min(1.0, rule["base_score"] * strength * (1.0 + 0.05 * bonus_hit))

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
    depths : shape (Z, pixels), band depth = 1 - R / R_continuum, >= 0;
        NaN where no continuum could be measured (the band has no neighbour
        inside its own coverage segment, or every such neighbour is nodata).
        NaN never satisfies a ``>=`` threshold test, so an unmeasurable band
        counts as "not detected" rather than as a spurious detection.

    Notes
    -----
    The continuum is taken only from bands in the same contiguous coverage
    segment as the feature (see :func:`_spectral_segments`), so it is never
    read from the far side of a spectral gap.
    """
    Z = len(wavelengths)
    depths = np.full_like(cube, np.nan)
    segments = _spectral_segments(wavelengths)
    for i in range(Z):
        mask = (np.abs(wavelengths - wavelengths[i]) < margin) & (segments == segments[i])
        mask[i] = False
        if not mask.any():
            # Isolated band: no shoulder, so no depth. Leaving this NaN is
            # the point -- substituting a continuum of 1.0 used to report
            # 1 - R here, which on a dark surface is a deep "absorption"
            # that does not exist.
            continue
        # All-NaN columns (true no-data pixels) are expected here and stay
        # NaN below; suppress the resulting "All-NaN slice" RuntimeWarning,
        # which np.errstate does not catch (it's raised via warnings, not a
        # fp exception).
        with np.errstate(all='ignore'), warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            local_max = np.nanmax(cube[mask], axis=0)
        with np.errstate(all='ignore'):
            usable = np.isfinite(local_max) & (local_max > 1e-10)
            depths[i] = np.where(
                usable, np.maximum(0.0, 1.0 - cube[i] / np.where(usable, local_max, 1.0)),
                np.nan)
    return depths


# Species that gate at least one rule's required_species -- these are the
# only ones worth computing scene-wide percentile statistics for (bonus/
# absent-only species are secondary signals, left on the absolute floor).
_REQUIRED_SPECIES_NAMES: set[str] = {
    sp for r in COMPOSITE_RULES for sp in r["required_species"]
}


def _accumulate_species_best_depth(depths_chunk: np.ndarray,
                                    sp_index: dict[str, list[int]],
                                    best_depth: dict[str, np.ndarray],
                                    pix_start: int, pix_end: int) -> None:
    """Fill best_depth[sp][pix_start:pix_end] with this chunk's per-pixel max
    depth within each required species' band set (NaN where that species has
    no bands in this cube)."""
    for sp, arr in best_depth.items():
        bands = sp_index.get(sp, [])
        if not bands:
            continue
        band_idxs = [i for i, _floor in bands]
        # All-NaN columns (true no-data pixels) are expected and yield NaN
        # here, same benign case as _compute_depths_numpy's nanmax guard.
        with np.errstate(invalid='ignore'), warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            arr[pix_start:pix_end] = np.nanmax(depths_chunk[band_idxs], axis=0)


# Validity codes written to the optional validity= raster. They separate the
# three things the class map alone conflates: a class was found, no class
# could be supported although the pixel was readable, and the pixel was not
# readable in the first place.
VALIDITY_NODATA = 0
VALIDITY_CLASSIFIED = 1
VALIDITY_UNCLASSIFIED = 2
VALIDITY_LOW_SIGNAL = 3
VALIDITY_SATURATED = 4
VALIDITY_LABELS = {
    VALIDITY_CLASSIFIED: "classified",
    VALIDITY_UNCLASSIFIED: "unclassified",
    VALIDITY_LOW_SIGNAL: "low signal",
    VALIDITY_SATURATED: "saturated",
}


def _detect_reflectance_scale(raster3d: str) -> float:
    """Divisor that brings the cube's stored values into 0-1 reflectance.

    Band depth is a ratio, so the engine's feature detection is indifferent to
    the storage scale; the absolute thresholds are not. Surface-reflectance
    products commonly ship scaled by 1000 or 10000 (the C2RD Pourrieres
    PARGE/ATCOR mosaic has a true range of -975 to 14981), and comparing such
    values against a threshold of 0.02 or 1.2 is meaningless. The scale is
    read from the cube's own range and reported, so a wrong guess is visible
    rather than silent.
    """
    from grass.exceptions import CalledModuleError

    try:
        info = gs.raster3d_info(raster3d)
        top = float(info["max"])
    except (CalledModuleError, KeyError, TypeError, ValueError):
        gs.warning("Could not read the cube's value range; assuming values "
                   "are already 0-1 reflectance.")
        return 1.0
    if not np.isfinite(top):
        return 1.0
    for limit, scale in ((1.5, 1.0), (150.0, 100.0), (1500.0, 1000.0)):
        if top <= limit:
            return scale
    return 10000.0


def _species_measurable(cube: np.ndarray,
                        sp_index: dict[str, list[tuple[int, float]]],
                        min_reflectance: float) -> dict[str, np.ndarray]:
    """Per species, whether this pixel is bright enough for a depth to mean
    something.

    Band depth is a ratio, ``1 - R / R_continuum``, so it says nothing about
    the absolute size of the dip. On a surface reflecting 2 % in the SWIR a
    depth of 0.3 is a reflectance difference of 0.006, at or below the noise
    of an atmospherically corrected imaging spectrometer. Measured on the
    C2RD Pourrieres targets, every apparent polyamide detection (N-H 1535 nm
    and amide 2049 nm depths of 0.12-0.31) came from a surface below 0.07
    SWIR reflectance, and the apparent depth correlated monotonically with
    darkness rather than with material. A species whose own bands are all
    darker than *min_reflectance* is therefore *indeterminate* here: neither
    detected nor absent, and no rule that needs it may fire.

    Returns one boolean array per species; True means measurable. A species
    with no bands in this sensor is absent from the result, which the caller
    treats as measurable so that the separate coverage test
    (:func:`unsupported_rules`) remains the single place that decides it.
    """
    out: dict[str, np.ndarray] = {}
    if min_reflectance <= 0:
        return out
    for sp, bands in sp_index.items():
        idxs = [i for i, _floor in bands]
        if not idxs:
            continue
        with np.errstate(invalid='ignore'), warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            brightest = np.nanmax(cube[idxs], axis=0)
        out[sp] = np.isfinite(brightest) & (brightest >= min_reflectance)
    return out


def _species_percentile_thresholds(best_depth: dict[str, np.ndarray],
                                    percentile: float) -> dict[str, float]:
    """Scene-wide percentile (0-100) of each species' best-depth distribution.

    percentile <= 0 disables the mechanism (returns an all-zero mapping,
    which combined via max() with the absolute floor is a no-op).
    """
    if percentile <= 0:
        return {sp: 0.0 for sp in best_depth}
    thresholds: dict[str, float] = {}
    for sp, arr in best_depth.items():
        valid = arr[~np.isnan(arr)]
        thresholds[sp] = float(np.percentile(valid, percentile)) if valid.size else 0.0
    return thresholds


def _build_band_species_index(
        wavelengths: np.ndarray) -> dict[str, list[tuple[int, float]]]:
    """For each DB species, the band indices in its range and their floors.

    A DB entry may carry its own ``min_depth``, which applies to the bands of
    that entry alone. A species whose bands differ in strength (the water OH
    overtones, for instance) needs this: one species-level threshold would be
    either too strict for the weak bands or too loose for the strong ones.
    Entries without ``min_depth`` contribute a floor of 0.0, so the rule's
    own override and the global minimum decide, exactly as before.
    """
    idx: dict[str, list[tuple[int, float]]] = {}
    for entry in SPECTRAL_FEATURES_DB:
        sp = entry["species"]
        lo, hi = entry["range"]
        floor = float(entry.get("min_depth", 0.0))
        hits = [(i, floor) for i, wl in enumerate(wavelengths) if lo <= wl <= hi]
        if hits:
            idx.setdefault(sp, []).extend(hits)
    # One entry per band, keeping the deepest floor where entries overlap.
    out: dict[str, list[tuple[int, float]]] = {}
    for sp, hits in idx.items():
        best: dict[int, float] = {}
        for i, floor in hits:
            best[i] = max(best.get(i, 0.0), floor)
        out[sp] = sorted(best.items())
    return out


# A band window is only "sampled" by a sensor when a band falls near its
# nominal centre, not merely somewhere inside the window. Without this test a
# sensor whose bands stop at 1482 nm counts as covering the 1400-1500 nm water
# window, and the engine then reads a band depth off the edge of a spectral
# gap -- which is how the 1940 nm water veto came to be evaluated on a cube
# that has no bands between 1784 and 1984 nm.
def _entry_is_sampled(wavelengths: np.ndarray, entry: dict) -> bool:
    lo, hi = entry["range"]
    tol = max(15.0, 0.4 * (hi - lo) / 2.0)
    centre = float(entry["center"])
    return bool(np.any(np.abs(wavelengths - centre) <= tol))


def _species_sampled(wavelengths: np.ndarray) -> dict[str, bool]:
    """Whether the sensor samples at least one band of each DB species."""
    sampled: dict[str, bool] = {}
    for entry in SPECTRAL_FEATURES_DB:
        sp = entry["species"]
        sampled[sp] = sampled.get(sp, False) or _entry_is_sampled(wavelengths, entry)
    return sampled


def unsupported_rules(wavelengths: np.ndarray) -> dict[str, str]:
    """Rules this sensor cannot decide, with the reason for each.

    A rule is unsupported when the sensor samples none of the bands of one of
    its required species, or none of the bands of one of its required-absent
    species. The second case is the dangerous one: an exclusion that cannot
    be tested would otherwise pass everywhere, turning a veto into a no-op
    and letting the rule fire on exactly the material it was meant to reject.
    Such a rule is disabled rather than silently relaxed. Bonus species are
    optional and never disable a rule.
    """
    sampled = _species_sampled(wavelengths)
    out: dict[str, str] = {}
    for rule in COMPOSITE_RULES:
        for sp in rule["required_species"]:
            if not sampled.get(sp, False):
                out[rule["name"]] = (
                    f"required species <{sp}> is outside this sensor's bands"
                )
                break
        else:
            for sp in rule.get("required_absent", []):
                if not sampled.get(sp, False):
                    out[rule["name"]] = (
                        f"exclusion species <{sp}> is outside this sensor's "
                        "bands, so the rule cannot be vetoed"
                    )
                    break
    return out


def warn_unsupported_rules(wavelengths: np.ndarray) -> dict[str, str]:
    """Report the unsupported rules to the user and return them."""
    unsupported = unsupported_rules(wavelengths)
    for name, reason in sorted(unsupported.items()):
        gs.warning(f"Rule <{name}> is disabled for this sensor: {reason}.")
    return unsupported


def score_composites_numpy(depths: np.ndarray, wavelengths: np.ndarray,
                            min_depth: float = 0.02,
                            species_percentile_thresholds: Optional[dict[str, float]] = None,
                            species_measurable: Optional[dict[str, np.ndarray]] = None
                            ) -> np.ndarray:
    """Score all composite rules per pixel using vectorized numpy.

    Parameters
    ----------
    depths      : shape (Z, pixels), band depths from _compute_depths_numpy
    wavelengths : shape (Z,)
    min_depth   : minimum depth to count a species as detected
    species_percentile_thresholds : optional species -> depth-value mapping,
        typically a percentile of that species' depth distribution across
        the whole scene. Combined with min_depth/per-rule overrides by
        taking the max, so a species is only "detected" where its depth is
        BOTH above the absolute floor AND an outlier relative to this
        scene's own baseline for that species -- guards against a broad,
        gently-sloped (not flat) continuum common to a material across the
        whole scene registering as a "feature" everywhere purely from local
        band-depth geometry, not a genuine localized absorption.

    species_measurable : optional species -> bool array, from
        :func:`_species_measurable`. False marks a species as *indeterminate*
        at that pixel: too dark for a band depth to carry information. A rule
        is not scored where any of its required or required-absent species is
        indeterminate, because "absent" and "unreadable" are different
        things and treating the second as the first is what lets a rule fire
        on a surface that cannot support it. Species missing from the mapping
        count as measurable.

    Returns
    -------
    scores : shape (n_rules, pixels), score 0–1 per rule per pixel. The value
        is the rule's base score scaled by how far its evidence sits above
        the detection threshold, so it reflects the strength of the match and
        not only which rule matched.
    """
    n_pixels = depths.shape[1]
    sp_index = _build_band_species_index(wavelengths)
    unsupported = unsupported_rules(wavelengths)
    n_rules = len(COMPOSITE_RULES)
    scores = np.zeros((n_rules, n_pixels), dtype=np.float32)
    pct_thresholds = species_percentile_thresholds or {}

    def _effective_threshold(sp: str, overrides: dict[str, float]) -> float:
        return max(overrides.get(sp, min_depth), pct_thresholds.get(sp, 0.0))

    # Pre-compute detection map: species → bool array (n_pixels,)
    def _hits(sp: str, threshold: float) -> np.ndarray:
        """Per-band detection, each band against the deeper of the two floors."""
        bands = sp_index.get(sp, [])
        if not bands:
            return np.zeros((0, n_pixels), dtype=bool)
        return np.stack([
            depths[i] >= max(threshold, floor) for i, floor in bands
        ])

    def _species_detected(sp: str, threshold: float) -> np.ndarray:
        hits = _hits(sp, threshold)
        if not hits.size:
            return np.zeros(n_pixels, dtype=bool)
        return np.any(hits, axis=0)

    def _species_n_detections(sp: str, threshold: float) -> np.ndarray:
        hits = _hits(sp, threshold)
        if not hits.size:
            return np.zeros(n_pixels, dtype=np.int32)
        return np.sum(hits, axis=0).astype(np.int32)

    measurable = species_measurable or {}

    def _species_measurable_here(sp: str) -> Optional[np.ndarray]:
        """None when nothing is known about this species' measurability."""
        return measurable.get(sp)

    def _species_evidence(sp: str, threshold: float) -> np.ndarray:
        """Strength of this species' evidence, in [0, 1].

        Two independent things make a detection convincing, and a band depth
        on its own reports neither: how far the deepest band sits above the
        threshold, and how much of the species' window agrees. A single band
        scraping over the threshold while its seven neighbours stay flat is a
        spike; the same depth across six of eight bands is an absorption
        band. They are combined as a geometric mean, so a high value needs
        both.
        """
        bands = sp_index.get(sp, [])
        if not bands:
            return np.zeros(n_pixels, dtype=np.float32)
        best = np.full(n_pixels, -np.inf, dtype=np.float32)
        detected = np.zeros(n_pixels, dtype=np.float32)
        for i, floor in bands:
            thresh = max(threshold, floor)
            with np.errstate(invalid='ignore'):
                over = np.where(np.isfinite(depths[i]), depths[i] - thresh, -np.inf)
            best = np.maximum(best, over.astype(np.float32))
            detected += (over >= 0.0).astype(np.float32)
        margin = np.clip(best / _EVIDENCE_DEPTH_SCALE, 0.0, 1.0)
        support = detected / float(len(bands))
        return np.sqrt(margin * support).astype(np.float32)

    for r_idx, rule in enumerate(COMPOSITE_RULES):
        if rule["name"] in unsupported:
            continue
        req = rule["required_species"]
        # Some rules need a stricter depth than the global default for a
        # specific species, when its generic detection window coincidentally
        # overlaps a common, unrelated material's ordinary spectral shape
        # (see acid_mine_drainage: a shallow jarosite-window dip is common in
        # ordinary hematite/goethite, not just genuine jarosite).
        overrides = rule.get("min_depth_overrides", {})

        # Count how many times each species is required
        req_counts: dict[str, int] = {}
        for sp in req:
            req_counts[sp] = req_counts.get(sp, 0) + 1

        # All required species must be detected sufficient times, and must be
        # measurable here at all: an indeterminate input decides nothing.
        all_met = np.ones(n_pixels, dtype=bool)
        margins: list[np.ndarray] = []
        for sp, cnt in req_counts.items():
            thresh = _effective_threshold(sp, overrides)
            if cnt == 1:
                all_met &= _species_detected(sp, thresh)
            else:
                all_met &= (_species_n_detections(sp, thresh) >= cnt)
            ok = _species_measurable_here(sp)
            if ok is not None:
                all_met &= ok
            margins.append(_species_evidence(sp, thresh))

        # Required absent species. An exclusion that cannot be read is not a
        # pass: the rule is withheld there, exactly as it is withheld when
        # the sensor does not sample the exclusion's bands at all.
        for sp in rule.get("required_absent", []):
            all_met &= ~_species_detected(sp, _effective_threshold(sp, overrides))
            ok = _species_measurable_here(sp)
            if ok is not None:
                all_met &= ok

        if not np.any(all_met):
            continue

        # Bonus species
        bonus_count = np.zeros(n_pixels, dtype=np.float32)
        for sp in rule.get("bonus_species", []):
            bonus_count += _species_detected(sp, _effective_threshold(sp, overrides)).astype(np.float32)

        # Calibrated confidence: the rule's base score is the ceiling, and
        # the depth margin of its weakest required species decides how much
        # of that ceiling this pixel earns. The geometric mean is used so one
        # species sitting at its threshold pulls the whole match down without
        # zeroing it, and so adding a required species cannot raise the
        # score. A match that only just clears every threshold keeps
        # _EVIDENCE_FLOOR of the base score; one well clear of all of them
        # keeps all of it, before the bonus factor.
        if margins:
            evidence = np.exp(np.mean(np.log(np.maximum(
                np.stack(margins), 1e-6)), axis=0)).astype(np.float32)
        else:
            evidence = np.ones(n_pixels, dtype=np.float32)
        strength = _EVIDENCE_FLOOR + (1.0 - _EVIDENCE_FLOOR) * evidence
        bonus_factor = 1.0 + 0.05 * bonus_count
        score = np.minimum(
            1.0, rule["base_score"] * strength * bonus_factor).astype(np.float32)
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
    validity_map = options.get('validity') or ''
    min_depth    = float(options['min_band_depth'])
    min_refl     = float(options.get('min_reflectance', '0.02') or '0.02')
    max_refl     = float(options.get('max_reflectance', '1.2') or '1.2')
    refl_scale   = float(options.get('reflectance_scale', '0') or '0')
    ndre_thresh  = float(options.get('ndre_threshold', '0.15') or '0.15')
    sand_albedo_thresh = float(options.get('sand_albedo_threshold', '0.25') or '0.25')
    species_percentile = float(options.get('species_percentile', '90') or '90')
    coordinates  = options.get('coordinates', '')
    only_valid   = flags['n']
    info_mode    = flags['i']
    point_mode   = flags['p']
    verbose      = flags['v']

    # Collect band metadata
    bands = get_band_info(input_map, only_valid=only_valid)

    # Absolute thresholds are expressed in 0-1 reflectance; the cube may be
    # stored scaled. Resolve the scale once and convert the thresholds into
    # the cube's own units, so every comparison below is in stored values.
    if refl_scale <= 0:
        refl_scale = _detect_reflectance_scale(input_map)
        gs.message(f"Reflectance scale detected: {refl_scale:g} "
                   f"(stored value / {refl_scale:g} = reflectance)")
    else:
        gs.verbose(f"Reflectance scale given: {refl_scale:g}")
    min_refl_stored = min_refl * refl_scale
    max_refl_stored = max_refl * refl_scale
    sand_albedo_stored = sand_albedo_thresh * refl_scale

    wavelengths = np.array([b['wavelength'] for b in bands], dtype=np.float64)
    fwhms       = np.array([b['fwhm']       for b in bands], dtype=np.float64)

    red_edge_idx = _nearest_band_index(wavelengths, _NDRE_RED_EDGE_NM)
    nir_idx = _nearest_band_index(wavelengths, _NDRE_NIR_NM)
    gs.verbose(
        f"NDRE pre-screen: red-edge band {red_edge_idx} "
        f"({wavelengths[red_edge_idx]:.1f}nm), NIR band {nir_idx} "
        f"({wavelengths[nir_idx]:.1f}nm), threshold={ndre_thresh}"
    )

    bright_lo, bright_hi = _SAND_BRIGHTNESS_RANGE_NM
    bright_mask = (wavelengths >= bright_lo) & (wavelengths <= bright_hi)
    gs.verbose(
        f"Quartz/feldspar sand fallback: broadband "
        f"{bright_lo:.0f}-{bright_hi:.0f}nm ({int(bright_mask.sum())} bands), "
        f"albedo threshold={sand_albedo_thresh}"
    )

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
        column = _spectrum_via_r3what(input_map, east, north)
        if column is not None:
            refl_values = [column[b['band'] - 1] for b in bands]
        else:
            gs.verbose("Falling back to per-band 2D slice extraction…")
            refl_values = []
            for z_idx, b in enumerate(bands):
                gs.percent(z_idx, len(bands), 2)
                tmp = extract_band(input_map, b['band'])
                val_str = gs.read_command(
                    'r.what', map=tmp, coordinates=f"{east},{north}", quiet=True
                ).split('|')[-1].strip()
                try:
                    refl_values.append(float(val_str))
                except ValueError:
                    refl_values.append(float('nan'))
            gs.percent(len(bands), len(bands), 2)

        spectrum: list[tuple[float, float, float]] = []
        band_refl: dict[int, float] = {}
        for z_idx, (b, refl) in enumerate(zip(bands, refl_values)):
            spectrum.append((b['wavelength'], b['fwhm'], refl))
            band_refl[z_idx] = refl

        ndre = _safe_ndre(np.array([band_refl[nir_idx]]),
                          np.array([band_refl[red_edge_idx]]))[0]
        gs.message(f"NDRE = {ndre:.3f} (threshold {ndre_thresh})")

        # Filter out nodata
        spectrum = [(c, fw, r) for c, fw, r in spectrum if not np.isnan(r)]
        if not spectrum:
            gs.fatal("All band values are nodata at the specified coordinates.")

        result = interpret_spectrum(
            spectrum, min_band_depth=min_depth, ndre=ndre,
            ndre_threshold=ndre_thresh,
            wavelengths=np.array([c for c, _fw, _r in spectrum]),
            min_reflectance=min_refl_stored)
        print_interpretation(result, verbose=verbose)
        return 0

    # ── Full-map mode ───────────────────────────────────────────────────────
    if not output_map and not conf_map and not validity_map:
        gs.fatal("Specify output=, confidence= and/or validity= for "
                 "full-map mode.")

    # Align the 3D computational region's depth window to the raster's own
    # native depth stack (bottom/top/tbres -- the caller's region may have
    # depths=1 by default) while preserving the caller's horizontal
    # extent/resolution. GRASS keeps a *separate* 3D region (rows3/cols3/
    # nsres3/ewres3) alongside the 2D one, and it is sticky across the
    # whole mapset (persisted in WIND, outliving this process). Never call
    # `g.region raster3d=<map>` -- it resets BOTH the 2D and 3D region to
    # the raster's full native footprint, silently discarding any smaller
    # area the caller had selected; merely restoring n/s/e/w/rows/cols
    # afterwards only fixes the 2D side, leaving 2D/3D mismatched. Instead
    # sync the 3D horizontal resolution to the already-correct 2D one
    # directly, and set the depth window from the raster's own metadata.
    depth_info = gs.parse_command('r3.info', map=input_map, flags='g')
    region = gs.region()
    gs.run_command('g.region',
                   nsres3=region['nsres'], ewres3=region['ewres'],
                   b=depth_info['bottom'], t=depth_info['top'],
                   tbres=depth_info['tbres'])

    nrows = int(region['rows'])
    ncols = int(region['cols'])
    n_pixels = nrows * ncols
    Z = len(bands)

    gs.message(f"Extracting {Z} bands ({nrows}×{ncols} pixels) via fast Z-slice…")
    tmp_names = []
    for z_idx, b in enumerate(bands):
        gs.percent(z_idx, Z, 2)
        tmp_names.append(extract_band(input_map, b['band']))
    gs.percent(Z, Z, 2)

    # Process the cube in row-chunks of ~TARGET_CHUNK_PIXELS pixels so peak
    # memory stays roughly constant regardless of scene size, instead of ever
    # materializing a full (bands, pixels) cube plus a same-size depths array
    # simultaneously (which for a full EMIT scene can require tens of GB).
    warn_unsupported_rules(wavelengths)

    gs.message("Computing band depths and scoring composites (tiled)…")
    TARGET_CHUNK_PIXELS = 300_000
    chunk_rows = max(1, TARGET_CHUNK_PIXELS // ncols)

    def _signal_quality(chunk: np.ndarray) -> tuple[np.ndarray, np.ndarray,
                                                     np.ndarray]:
        """Classify each pixel's readability: nodata, saturated, low signal.

        A pixel whose whole spectrum is darker than min_reflectance carries
        no usable absorption information at all, and one with a band above
        max_reflectance has been clipped by the sensor or driven out of range
        by the atmospheric correction (the white glass-fibre panel on the
        C2RD scene reads 113 % at 550 nm and nodata through the SWIR). Both
        are reported in the validity map and left unclassified, instead of
        being passed off as an ordinary no-match.
        """
        nodata = np.all(np.isnan(chunk), axis=0)
        with np.errstate(invalid='ignore'), warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            brightest = np.nanmax(chunk, axis=0)
        saturated = np.isfinite(brightest) & (brightest > max_refl_stored) & ~nodata
        low = (
            np.isfinite(brightest) & (brightest < min_refl_stored) & ~nodata & ~saturated
            if min_refl_stored > 0 else np.zeros_like(nodata)
        )
        return nodata, saturated, low

    class_map = np.zeros(n_pixels, dtype=np.int32)
    conf_arr = np.zeros(n_pixels, dtype=np.float32)
    nodata_mask = np.zeros(n_pixels, dtype=bool)
    validity_arr = np.zeros(n_pixels, dtype=np.int32)
    n_ndre_screened = 0
    n_sand_fallback = 0
    n_low_signal = 0
    n_saturated = 0

    try:
        from grass.pygrass.raster import RasterRow
        handles = [RasterRow(name) for name in tmp_names]
        for h in handles:
            h.open()
        use_pygrass = True
    except Exception:
        gs.warning("grass.pygrass unavailable — falling back to loading the "
                   "full cube at once (higher memory use).")
        use_pygrass = False

    if use_pygrass:
        try:
            n_chunks = (nrows + chunk_rows - 1) // chunk_rows

            def _read_chunk(row_start: int, row_end: int) -> np.ndarray:
                n_pix = (row_end - row_start) * ncols
                chunk = np.empty((Z, n_pix), dtype=np.float32)
                for z_idx, h in enumerate(handles):
                    for r in range(row_start, row_end):
                        row_arr = np.array(h[r], dtype=np.float32)
                        offset = (r - row_start) * ncols
                        chunk[z_idx, offset:offset + ncols] = row_arr
                return chunk

            # Pass 1: accumulate each required species' best (max) depth per
            # pixel across the whole scene, to derive scene-relative
            # percentile thresholds (see score_composites_numpy docstring)
            # before any classification decision is made.
            sp_index = _build_band_species_index(wavelengths)
            best_depth = {sp: np.full(n_pixels, np.nan, dtype=np.float32)
                          for sp in _REQUIRED_SPECIES_NAMES}
            if species_percentile > 0:
                gs.message(f"Pass 1/2: scene depth statistics "
                           f"(percentile={species_percentile})…")
                for c_idx, row_start in enumerate(range(0, nrows, chunk_rows)):
                    gs.percent(c_idx, n_chunks, 2)
                    row_end = min(row_start + chunk_rows, nrows)
                    chunk = _read_chunk(row_start, row_end)
                    depths_chunk = _compute_depths_numpy(chunk, wavelengths,
                                                          margin=_CONTINUUM_MARGIN)
                    pix_start = row_start * ncols
                    pix_end = row_end * ncols
                    _accumulate_species_best_depth(depths_chunk, sp_index,
                                                    best_depth, pix_start, pix_end)
                    del chunk, depths_chunk
                gs.percent(n_chunks, n_chunks, 2)
            species_pct_thresholds = _species_percentile_thresholds(
                best_depth, species_percentile)
            del best_depth
            if species_percentile > 0:
                for sp, thresh in sorted(species_pct_thresholds.items(),
                                         key=lambda kv: -kv[1])[:10]:
                    gs.verbose(f"  species_percentile[{sp}] = {thresh:.4f}")

            # Pass 2: classify, gating required-species detection by both the
            # absolute floor and the scene-relative percentile threshold.
            gs.message("Pass 2/2: scoring composites…")
            for c_idx, row_start in enumerate(range(0, nrows, chunk_rows)):
                gs.percent(c_idx, n_chunks, 2)
                row_end = min(row_start + chunk_rows, nrows)
                n_pix_chunk = (row_end - row_start) * ncols
                chunk = _read_chunk(row_start, row_end)

                depths_chunk = _compute_depths_numpy(chunk, wavelengths,
                                                      margin=_CONTINUUM_MARGIN)
                scores_chunk = score_composites_numpy(
                    depths_chunk, wavelengths, min_depth=min_depth,
                    species_percentile_thresholds=species_pct_thresholds,
                    species_measurable=_species_measurable(
                        chunk, sp_index, min_refl_stored))

                # NDRE pre-screen: bar chlorophyll-dependent rules (vegetation,
                # algae blooms) on pixels that don't actually look vegetated,
                # regardless of what the generic band-depth matching found.
                ndre_chunk = _safe_ndre(chunk[nir_idx], chunk[red_edge_idx])
                veg_blocked = ndre_chunk < ndre_thresh
                if veg_blocked.any():
                    for gi in _CHLOROPHYLL_GATED_RULE_INDICES:
                        scores_chunk[gi, veg_blocked] = 0.0
                    n_ndre_screened += int(veg_blocked.sum())

                best_idx = np.argmax(scores_chunk, axis=0)
                best_score = scores_chunk[best_idx, np.arange(n_pix_chunk)]

                cls = (best_idx + 1).astype(np.int32)
                cls[best_score < 0.1] = 0
                conf = best_score.astype(np.float32)

                nodata_chunk, saturated_chunk, low_chunk = _signal_quality(chunk)
                cls[nodata_chunk | saturated_chunk | low_chunk] = 0
                conf[nodata_chunk] = np.nan
                conf[saturated_chunk | low_chunk] = 0.0
                n_saturated += int(saturated_chunk.sum())
                n_low_signal += int(low_chunk.sum())

                # Quartz/feldspar sand fallback: no diagnostic feature exists
                # for these, so pixels left unclassified (no rule fired),
                # already screened as non-vegetated, and bright are assigned
                # here by elimination rather than a positive species match.
                broadband_chunk = _broadband_reflectance(chunk, bright_mask)
                sand_fallback = (
                    (cls == 0) & ~nodata_chunk & ~saturated_chunk & ~low_chunk
                    & (ndre_chunk < ndre_thresh)
                    & (broadband_chunk > sand_albedo_stored)
                )
                cls[sand_fallback] = QUARTZOFELDSPATHIC_SAND_CLASS_IDX
                conf[sand_fallback] = QUARTZOFELDSPATHIC_SAND_CONFIDENCE
                n_sand_fallback += int(sand_fallback.sum())

                pix_start = row_start * ncols
                pix_end = row_end * ncols
                validity_chunk = np.where(
                    nodata_chunk, VALIDITY_NODATA,
                    np.where(saturated_chunk, VALIDITY_SATURATED,
                             np.where(low_chunk, VALIDITY_LOW_SIGNAL,
                                      np.where(cls > 0, VALIDITY_CLASSIFIED,
                                               VALIDITY_UNCLASSIFIED))))

                class_map[pix_start:pix_end] = cls
                conf_arr[pix_start:pix_end] = conf
                nodata_mask[pix_start:pix_end] = nodata_chunk
                validity_arr[pix_start:pix_end] = validity_chunk

                del chunk, depths_chunk, scores_chunk
            gs.percent(n_chunks, n_chunks, 2)
        finally:
            for h in handles:
                h.close()
    else:
        # Fallback: original full-cube path (float32 instead of float64 to
        # at least halve peak memory versus the pre-fix behaviour).
        cube = np.full((Z, n_pixels), np.nan, dtype=np.float32)
        for z_idx, tmp in enumerate(tmp_names):
            ascii_out = gs.tempfile()
            gs.run_command('r.out.ascii', input=tmp, output=ascii_out,
                           null_value='nan', quiet=True)
            arr = np.loadtxt(ascii_out, dtype=np.float32).ravel()
            os.unlink(ascii_out)
            if arr.size == n_pixels:
                cube[z_idx] = arr
            else:
                gs.warning(f"Band {bands[z_idx]['band']}: unexpected array "
                           f"size {arr.size}, expected {n_pixels}")

        depths = _compute_depths_numpy(cube, wavelengths, margin=_CONTINUUM_MARGIN)

        sp_index = _build_band_species_index(wavelengths)
        best_depth = {sp: np.full(n_pixels, np.nan, dtype=np.float32)
                      for sp in _REQUIRED_SPECIES_NAMES}
        _accumulate_species_best_depth(depths, sp_index, best_depth, 0, n_pixels)
        species_pct_thresholds = _species_percentile_thresholds(
            best_depth, species_percentile)
        del best_depth

        scores = score_composites_numpy(
            depths, wavelengths, min_depth=min_depth,
            species_percentile_thresholds=species_pct_thresholds,
            species_measurable=_species_measurable(cube, sp_index, min_refl_stored))

        ndre_full = _safe_ndre(cube[nir_idx], cube[red_edge_idx])
        veg_blocked = ndre_full < ndre_thresh
        if veg_blocked.any():
            for gi in _CHLOROPHYLL_GATED_RULE_INDICES:
                scores[gi, veg_blocked] = 0.0
            n_ndre_screened += int(veg_blocked.sum())

        best_idx = np.argmax(scores, axis=0)
        best_score = scores[best_idx, np.arange(n_pixels)]
        class_map = (best_idx + 1).astype(np.int32)
        class_map[best_score < 0.1] = 0
        conf_arr = best_score.astype(np.float32)
        nodata_mask, saturated_full, low_full = _signal_quality(cube)
        class_map[nodata_mask | saturated_full | low_full] = 0
        conf_arr[nodata_mask] = np.nan
        conf_arr[saturated_full | low_full] = 0.0
        n_saturated += int(saturated_full.sum())
        n_low_signal += int(low_full.sum())

        broadband_full = _broadband_reflectance(cube, bright_mask)
        sand_fallback = (
            (class_map == 0) & ~nodata_mask & ~saturated_full & ~low_full
            & (ndre_full < ndre_thresh)
            & (broadband_full > sand_albedo_stored)
        )
        class_map[sand_fallback] = QUARTZOFELDSPATHIC_SAND_CLASS_IDX
        conf_arr[sand_fallback] = QUARTZOFELDSPATHIC_SAND_CONFIDENCE
        n_sand_fallback += int(sand_fallback.sum())

        validity_arr = np.where(
            nodata_mask, VALIDITY_NODATA,
            np.where(saturated_full, VALIDITY_SATURATED,
                     np.where(low_full, VALIDITY_LOW_SIGNAL,
                              np.where(class_map > 0, VALIDITY_CLASSIFIED,
                                       VALIDITY_UNCLASSIFIED)))).astype(np.int32)

        del cube, depths, scores

    class_map = class_map.reshape(nrows, ncols)
    conf_arr = conf_arr.reshape(nrows, ncols)
    nodata_mask = nodata_mask.reshape(nrows, ncols)
    validity_arr = validity_arr.reshape(nrows, ncols)

    if n_ndre_screened:
        gs.message(
            f"NDRE pre-screen (threshold={ndre_thresh}): barred "
            f"{n_ndre_screened} pixel(s) from chlorophyll-dependent classes "
            f"({', '.join(sorted(_CHLOROPHYLL_GATED_RULE_NAMES))})"
        )
    if n_low_signal:
        gs.message(
            f"Low signal (all bands below {min_refl:g} reflectance, "
            f"{min_refl_stored:g} stored): "
            f"{n_low_signal} pixel(s) left unclassified; a band depth there "
            f"would be a ratio of noise"
        )
    if n_saturated:
        gs.message(
            f"Saturated (a band above {max_refl:g} reflectance, "
            f"{max_refl_stored:g} stored): "
            f"{n_saturated} pixel(s) left unclassified"
        )
    if n_sand_fallback:
        gs.message(
            f"Quartz/feldspar sand fallback (albedo>{sand_albedo_thresh}): "
            f"assigned {n_sand_fallback} otherwise-unclassified pixel(s) to "
            f"'{QUARTZOFELDSPATHIC_SAND_LABEL}' (confidence="
            f"{QUARTZOFELDSPATHIC_SAND_CONFIDENCE}, elimination-based)"
        )

    # Write output rasters
    if _RAS3D:
        import ras3d as _r3
        import ras3d_write as _r3w
        _h = _r3.open_cube(input_map)
        if output_map:
            gs.message(f"Writing class map → {output_map}")
            _r3w.write_raster2d(_r3w.outpath(output_map), class_map.astype(np.float32), _h)
        if conf_map:
            gs.message(f"Writing confidence map → {conf_map}")
            _r3w.write_raster2d(_r3w.outpath(conf_map), conf_arr, _h)
        if validity_map:
            gs.message(f"Writing validity map → {validity_map}")
            _r3w.write_raster2d(_r3w.outpath(validity_map),
                                validity_arr.astype(np.float32), _h)
        _r3.close_cube(_h)
        return 0

    reg = gs.region()
    _ascii_header = (
        f"north: {reg['n']}\n"
        f"south: {reg['s']}\n"
        f"east: {reg['e']}\n"
        f"west: {reg['w']}\n"
        f"rows: {nrows}\n"
        f"cols: {ncols}\n"
    )

    def _write_raster(name, array, rtype, null_val):
        ascii_tmp = gs.tempfile()
        with open(ascii_tmp, 'w') as _fh:
            _fh.write(_ascii_header)
            if rtype == 'int':
                np.savetxt(_fh, array, fmt='%d', delimiter=' ')
            else:
                np.savetxt(_fh, array, fmt='%.4f', delimiter=' ')
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
                           f"0=unknown, 1..{len(COMPOSITE_RULES)}=composite rules, "
                           f"{QUARTZOFELDSPATHIC_SAND_CLASS_IDX}="
                           f"{QUARTZOFELDSPATHIC_SAND_LABEL} (elimination-based "
                           f"fallback, not a composite rule). Input: {input_map}"))

    if validity_map:
        gs.message(f"Writing validity map → {validity_map}")
        _write_raster(validity_map, validity_arr, 'int', VALIDITY_NODATA)
        cats = "\n".join(f"{code}:{label}"
                          for code, label in sorted(VALIDITY_LABELS.items()))
        gs.write_command('r.category', map=validity_map,
                         rules='-', separator=':', stdin=cats, quiet=True)
        gs.run_command('r.support', map=validity_map,
                       title="Spectral interpretation - per-pixel validity",
                       units="validity code",
                       description=(
                           f"i.hyper.spectroscopy validity map. "
                           f"{VALIDITY_CLASSIFIED}=classified, "
                           f"{VALIDITY_UNCLASSIFIED}=unclassified (readable, "
                           f"no rule fired), {VALIDITY_LOW_SIGNAL}=low signal "
                           f"(every band below {min_refl:g} reflectance), "
                           f"{VALIDITY_SATURATED}=saturated (a band above "
                           f"{max_refl:g}); NULL where the cube has no data. "
                           f"Input: {input_map}"))
        gs.raster_history(validity_map, overwrite=True)

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
