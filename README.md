# i.hyper.spectroscopy

**GRASS GIS module — Per-pixel spectral feature interpretation of hyperspectral
imagery (400–2500 nm VNIR+SWIR)**

> **Database:** 163 absorption features · 58 composite material classes · 51 automated tests

Part of the [i.hyper](../README.md) module family for hyperspectral data
processing in GRASS GIS.

---

## Overview

`i.hyper.spectroscopy` interprets the reflectance spectrum at every pixel of a
hyperspectral 3D raster and assigns ranked material hypotheses based on detected
absorption features.

Rather than mapping spectral indices to fixed class thresholds, the engine works
in three layers:

1. **Feature layer** — identify absorption features from local band depths
   (1 − R_band / R_continuum), with a ±150 nm sliding continuum window
2. **Chromophore layer** — match each feature to candidate absorbing species
   (chlorophyll, Fe³⁺, Al-OH, CO₃²⁻, CH overtones, …) from a **163-entry**
   physics database spanning 400–2500 nm
3. **Evidence layer** — score **58 composite material hypotheses** by how many
   independent species co-occur (e.g., chlorophyll + red-edge → green vegetation;
   Al-OH + CO₃²⁻ → calcite/dolomite)

Two outputs are produced per pixel:

| Output | Type | Content |
|--------|------|---------|
| `output` | integer CELL | Dominant class code (0 = unknown, 1–58 = named class) |
| `confidence` | float FCELL | Best-hypothesis score 0–1 |

---

## Class codes

| Code | Name | Category |
|------|------|----------|
| 0 | unknown | — |
| 1 | green\_vegetation | biologic |
| 2 | stressed\_vegetation | biologic |
| 3 | dry\_plant\_litter | biologic |
| 4 | oxyhemoglobin\_blood | biologic |
| 5 | deoxyhemoglobin\_blood | biologic |
| 6 | kaolinite | mineral |
| 7 | smectite\_montmorillonite | mineral |
| 8 | muscovite\_illite | mineral |
| 9 | calcite | mineral |
| 10 | dolomite | mineral |
| 11 | gypsum\_evaporite | mineral |
| 12 | liquid\_water | material |
| 13 | ice\_snow | geological |
| 14 | olivine | mineral |
| 15 | pyroxene\_orthopyroxene | mineral |
| 16 | hematite\_iron\_oxide | mineral |
| 17 | goethite | mineral |
| 18 | neodymium\_REE | mineral |
| 19 | crude\_oil\_petroleum | organic |
| 20 | coal | geological |
| 21 | polyethylene\_plastic | synthetic |
| 22 | protein\_rich\_material | biologic |
| 23 | cyanobacteria\_algae | biologic |
| 24 | serpentinite | mineral |
| 25 | advanced\_argillic\_alteration | geological |
| 26 | propylitic\_alteration | geological |
| 27 | chrysotile\_asbestos | mineral |
| 28 | starch\_food\_grain | biologic |
| 29 | anthocyanin\_stress | biologic |
| 30 | phycoerythrin\_algae | biologic |
| 31 | cyanobacteria\_bloom | biologic |
| 32 | diatom\_bloom | biologic |
| 33 | ferrihydrite | mineral |
| 34 | lepidocrocite\_soil | mineral |
| 35 | acid\_mine\_drainage | geological |
| 36 | siderite | mineral |
| 37 | REE\_erbium | mineral |
| 38 | REE\_dysprosium | mineral |
| 39 | REE\_ytterbium | mineral |
| 40 | REE\_holmium | mineral |
| 41 | pyrophyllite | mineral |
| 42 | gibbsite\_laterite | mineral |
| 43 | boehmite\_bauxite | mineral |
| 44 | halloysite | mineral |
| 45 | talc | mineral |
| 46 | epsomite | mineral |
| 47 | kieserite | mineral |
| 48 | szomolnokite | mineral |
| 49 | prehnite | mineral |
| 50 | schorl\_tourmaline | mineral |
| 51 | crocidolite\_asbestos | mineral |
| 52 | ankerite | mineral |
| 53 | PET\_plastic | synthetic |
| 54 | PVC\_plastic | synthetic |
| 55 | polystyrene | synthetic |
| 56 | rubber\_polyisoprene | synthetic |
| 57 | concrete\_urban | urban |
| 58 | asphalt\_bitumen | urban |

---

## Absorption feature database

The built-in database covers **163 entries** grouped by physical mechanism:

| Group | Species | Key wavelengths (nm) |
|-------|---------|----------------------|
| Vegetation pigments | Chlorophyll a+b, carotenoids, phycocyanin | 430, 480, 620, 665, 710 |
| Stress pigments | Anthocyanins, zeaxanthin/PRI | 531, 540 |
| Aquatic pigments | Phycoerythrin, allophycocyanin, fucoxanthin | 490, 495, 545, 650 |
| Hemoglobin | HbO₂ (Soret, Q-bands), Hb, MetHb | 415, 542, 560, 577, 630 |
| Iron oxides (common) | Fe³⁺ hematite/goethite, Fe²⁺ pyroxene/olivine | 490, 530, 700, 900, 1000, 2000 |
| Iron oxides (expanded) | Ferrihydrite, lepidocrocite, siderite FeCO₃ | 430, 475, 750, 1000, 1080, 2340 |
| Rare earth elements | Nd³⁺ (×4), Sm³⁺/Er³⁺, Pr³⁺ | 525, 580, 745, 800, 867, 940 |
| REE expanded | Dy³⁺, Yb³⁺, Er³⁺ (×3), Ho³⁺, Eu³⁺, Tb³⁺, Ce³⁺ | 400, 450, 490, 535, 537, 650, 910, 975, 1260, 1530 |
| Other metals | Cu²⁺, Cr³⁺, Co²⁺, Mn²⁺ | 435, 460, 550, 630, 700 |
| Water / ice | Liquid H₂O, structural OH, ice | 970, 1200, 1380, 1450, 1940, 2000 |
| Clay minerals (common) | Al-OH, Mg-OH, NH₄⁺ | 1410, 2160, 2200, 2250, 2310 |
| Clay minerals (expanded) | Pyrophyllite, gibbsite, boehmite, halloysite, phengite, palygorskite, sepiolite, talc | 1390, 2110, 2165, 2175, 2205, 2215, 2263, 2315, 2387, 2390 |
| Carbonates | CO₃²⁻ calcite/dolomite/ankerite/siderite | 2120, 2320, 2330, 2335, 2350, 2500 |
| Sulphates (common) | Gypsum, alunite, jarosite | 1450, 1490, 1750, 2210, 2265 |
| Sulphates (expanded) | Epsomite, bassanite, kieserite, szomolnokite | 1000, 1490, 1640, 1940, 2165, 2200, 2220 |
| Organics | Cellulose, lignin, lipid, protein, starch | 1210, 1480, 1680, 1720, 2054, 2100 |
| Petroleum / coal | CH overtones, aromatic CH | 1210, 1720, 2300, 2310 |
| Polymers (common) | Polyethylene, polypropylene, nylon | 1150, 1215, 1660, 2310 |
| Polymers (expanded) | PET, PVC, polystyrene, rubber | 1640, 1680, 1720, 1730, 2030, 2170, 2240 |
| Silicates / feldspars | Si-OH, Fe²⁺ feldspar | 1250, 1380, 2200 |
| Sheet silicates | Mica, amphibole, serpentine, chlorite, epidote | 2230, 2250, 2310, 2340, 2370 |
| Geological minerals | Prehnite, schorl/tourmaline, crocidolite | 440, 1490, 2120, 2243, 2320 |
| Urban materials | Portlandite (concrete) | 1460 |
| Atmospheric refs | O₂ A/B bands, H₂O vapour columns | 690, 760, 820, 940, 1140 |
| Soil nutrients | Nitrate NO₃⁻, ammonium NH₄⁺ (2nd band) | 1410, 1560, 2050 |

---

## Usage

```bash
# Full-map classification (dominant class + confidence)
i.hyper.spectroscopy input=enmap_reflectance \
                     output=enmap_class \
                     confidence=enmap_conf

# Only valid bands; lower detection threshold
i.hyper.spectroscopy input=enmap_reflectance \
                     output=enmap_class \
                     confidence=enmap_conf \
                     min_band_depth=0.01 \
                     -n

# Info mode: inspect band wavelengths without processing
i.hyper.spectroscopy input=enmap_reflectance -i

# Point mode: interpret spectrum at a single location
i.hyper.spectroscopy input=enmap_reflectance \
                     coordinates=500000,5400000 \
                     -p -v
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `input` | — | Input hyperspectral 3D raster |
| `output` | — | Output class map (integer, 0–28) |
| `confidence` | — | Output confidence map (float 0–1) |
| `min_band_depth` | 0.02 | Minimum band depth to count as detected absorption |
| `coordinates` | — | `east,north` for point-mode interpretation |

### Flags

| Flag | Effect |
|------|--------|
| `-n` | Use only bands marked `valid=1` in metadata |
| `-i` | Info mode — print band/wavelength table and exit |
| `-p` | Point mode — interpret spectrum at `coordinates=` |
| `-v` | Verbose — print full hypothesis ranking (point mode) |

---

## Typical workflow

```bash
# 1. Import hyperspectral data
i.hyper.import input=/data/ENMAP_L2A.hdf \
               product=enmap \
               output=enmap

# 2. Atmospheric correction
i.hyper.smac input=enmap output=enmap_boa

# 3. (Optional) Continuum removal for sharper feature depths
i.hyper.continuum input=enmap_boa output=enmap_cr

# 4. Check sensor coverage
i.hyper.spectroscopy input=enmap_cr -i

# 5. Point-mode exploration at a known mineral outcrop
i.hyper.spectroscopy input=enmap_cr \
                     coordinates=450320,5317880 \
                     -p -v

# 6. Full classification
i.hyper.spectroscopy input=enmap_cr \
                     output=enmap_class \
                     confidence=enmap_conf \
                     -n

# 7. Visualise
d.rast map=enmap_class
d.legend map=enmap_class
r.univar map=enmap_conf zones=enmap_class
```

---

## Spectral rules of thumb

```
FWHM < 10 nm    → electronic transition (f→f REE, O₂, sharp molecular)
FWHM 10–50 nm   → crystal-field (Fe, Cu, Cr) or chlorophyll
FWHM 50–150 nm  → vibrational overtone / combination (OH, CH, NH, CO₃)
2200 nm alone   → Al-OH clay (kaolinite / illite / smectite)
1450 + 1940 nm  → liquid H₂O present
670 + red-edge  → photosynthetic vegetation
Band I ~900 + Band II ~2000 → pyroxene
Narrow 580+745+800 nm → REE Nd³⁺
1720 + 2310 nm (narrow FWHM) → petroleum / hydrocarbon
```

---

## Testing

The pure-Python interpretation engine (no GRASS session required) is covered by
`test_suite.py`:

```bash
python -m pytest test_suite.py -v
```

51 tests across 9 classes:

| Class | What is tested |
|-------|---------------|
| `TestDatabaseIntegrity` | Required fields on all DB entries; wavelength range sanity; all `required_species`/`bonus_species` exist in the DB; `base_score` in [0,1]; `CLASS_LABELS` consistency |
| `TestBandDepth` | Normal dip, zero reflectance, no neighbours (continuum = 1.0), feature = continuum → 0, feature brighter than continuum → clamped to 0 |
| `TestInterpretSpectrumEdgeCases` | Empty input, out-of-range wavelengths, sub-threshold FWHM, flat spectrum, no 60 % confidence flag, reflectance clamping, asymmetry classification |
| `TestInterpretSpectrumMaterials` | Synthetic spectra for green vegetation, liquid water, kaolinite, stressed vegetation (chlorophyll without red-edge), atmospheric contamination flag, hypothesis sort order, 60 % dominance threshold |
| `TestComputeDepthsNumpy` | Output shape, flat → zero depths, single-dip exact depth, isolated bands (no neighbours), NaN handling, all depths ≥ 0 |
| `TestBuildBandSpeciesIndex` | Known species at correct band indices, out-of-range wavelengths → empty, deduplication/sort |
| `TestScoreCompositesNumpy` | Shape, vegetation pixel ≥ 0.60 / flat pixel = 0, all scores in [0,1], `required_absent` exclusion |
| `TestChromophoreScoring` | Deeper feature → higher confidence; multi-band Nd³⁺ → higher confidence than single band |
| `TestNewMaterialClasses` | Pyrophyllite vs kaolinite disambiguation; ferrihydrite doublet; PVC/PET plastic; gibbsite doublet; crocidolite vs chrysotile; acid mine drainage; talc doublet; phycoerythrin doublet; DB integrity for all 30 new rules |

The test file mocks `grass.script` and loads `i.hyper.spectroscopy.py` via
`importlib`, so no GRASS installation is needed.

---

## Extending the engine

### Add an absorption feature to the database

```python
SPECTRAL_FEATURES_DB.append({
    "center": 2390,
    "range": (2375, 2405),
    "species": "Mg-OH_brucite",
    "mechanism": "Mg-OH_combination",
    "tags": ["mineral", "geological"],
    "notes": "Brucite Mg(OH)2; ultramafic alteration"
})
```

### Add a composite rule

```python
COMPOSITE_RULES.append({
    "name": "brucite",
    "category": "mineral",
    "required_species": ["Mg-OH_brucite", "Mg-OH_Fe-OH"],
    "base_score": 0.82,
    "notes": "Brucite: 2390 + 2250 Mg-OH"
})
```

---

## Sensor compatibility

| Sensor | VNIR | SWIR | Notes |
|--------|:----:|:----:|-------|
| HySpex Mjolnir VS-620 (400–2500 nm) | ✓ | ✓ | Full 28-class assessment |
| EnMAP (420–2450 nm, 224 bands) | ✓ | ✓ | Full assessment |
| PRISMA (400–2500 nm, 238 bands) | ✓ | ✓ | Full assessment |
| EMIT (380–2500 nm, 285 bands) | ✓ | ✓ | Full assessment |
| Sentinel-2 MSI (12 bands) | ✓ | partial | Vegetation + broad Fe³⁺ only |
| WV3-SWIR (8 bands, 1195–2365 nm) | ✗ | ✓ | Mineral classes only |

Sensors without coverage for a spectral region score zero for rules requiring
bands in that region; the remaining rules still apply.

---

## Disambiguation tips

| FWHM | Likely mechanism | Examples |
|------|-----------------|----------|
| < 10 nm | f→f electronic, O₂ | REE (Nd³⁺), atmospheric O₂ |
| 10–50 nm | Crystal-field, π→π* | Fe oxides, chlorophyll, Cu²⁺ |
| 50–150 nm | Vibrational overtone | OH, CH, NH, CO₃ |
| > 150 nm | Broad mixed-phase | Aqueous solution, amorphous |

Critical overlap zones:

| Zone (nm) | Competing species | Disambiguate by |
|-----------|------------------|-----------------|
| 2280–2380 | CO₃, CH hydrocarbon, Mg-OH | FWHM; rock vs. organic context |
| 1380–1420 | Structural OH, Al-OH, aromatic CH | FWHM; paired features |
| 1440–1500 | Gypsum H₂O, cellulose OH, liquid water | 1940 nm present → water; 1750 → gypsum |
| 700–730 | Fe³⁺ CT, red-edge, Cu²⁺ | Presence of 530 (Fe) or 665 (Chl) |

---

## References

- Hunt, G. R. (1977). Spectral signatures of particulate minerals in the
  visible and near infrared. *Geophysics*, 42(3), 501–513.
- Clark, R. N., et al. (1990). High spectral resolution reflectance
  spectroscopy of minerals. *Journal of Geophysical Research*, 95(B8),
  12653–12680.
- Kokaly, R. F., et al. (2017). USGS Spectral Library Version 7.
  *USGS Data Series* 1035. <https://doi.org/10.3133/ds1035>
- Clevers, J. G. P. W., & Kooistra, L. (2012). Using hyperspectral remote
  sensing data for retrieving canopy chlorophyll and nitrogen content.
  *IEEE JSTARS*, 5(2), 574–583.
- Vane, G., & Goetz, A. F. H. (1993). Terrestrial imaging spectrometry:
  Current status, future trends. *Remote Sensing of Environment*, 44(2–3),
  117–126.

## See also

- [i.hyper.geology](../i.hyper.geology/README.md) — geological classification
- [i.hyper.continuum](../i.hyper.continuum/README.md) — continuum removal
- [i.hyper.indices](../i.hyper.indices/README.md) — spectral indices (86+)
- [i.hyper.atcorr](../i.hyper.atcorr/README.md) — atmospheric correction (6SV2.1)
- [i.hyper.smac](../i.hyper.smac/README.md) — atmospheric correction (SMAC)
- GRASS manual page: [i.hyper.spectroscopy.html](i.hyper.spectroscopy.html)

---

## License

Unlicense — see [LICENSE](LICENSE)

## Authors

Created for the i.hyper module family (2025).
