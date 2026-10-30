## DESCRIPTION

*i.hyper.spectroscopy* performs per-pixel spectral feature interpretation
of hyperspectral imagery imported as 3D raster maps (`raster_3d`) by
[i.hyper.import](i.hyper.import.html).

The module reads wavelength metadata from hyperspectral 3D raster bands,
computes local band depths across the VNIR-SWIR range (400–2500 nm), matches
detected absorption features to a physics-based database of **163 absorbing
species entries**, and scores **58 composite material hypotheses** per pixel.
Two output rasters are produced: a dominant-class integer map and a confidence map.

*i.hyper.spectroscopy* is part of the **i.hyper** module family.

### Interpretation pipeline

The engine works in three layers:

1. **Feature layer** — For each band, the local band depth is estimated as
   `max(0, 1 − R_band / R_continuum)`, where the continuum is the maximum
   reflectance within a ±150 nm sliding window.  Bands with depth below
   `min_band_depth` are discarded.
2. **Chromophore layer** — Detected features are matched against the built-in
   absorption feature database.  Each match scores the candidate absorbing
   species by band-depth strength and detection multiplicity.
3. **Evidence layer** — Composite signatures (rules requiring specific
   combinations of species) are evaluated per pixel.  Rules that pass all
   required-species tests receive a score equal to their `base_score` plus a
   bonus for each optional bonus species that is also detected.

The class code with the highest score is written to the output map.  Pixels
where all composite scores fall below 0.10 receive class 0 (unknown).

### Class codes

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

### Absorption feature database

The built-in database covers **163 entries** across the full VNIR-SWIR range:

| Group | Key species | Diagnostic wavelengths (nm) |
|-------|-------------|------------------------------|
| Vegetation pigments | Chlorophyll a+b, carotenoids, phycocyanin | 430, 480, 620, 665, 710 |
| Stress pigments | Anthocyanins, zeaxanthin/PRI | 531, 540 |
| Aquatic pigments | Phycoerythrin, allophycocyanin, fucoxanthin | 490, 495, 545, 650 |
| Hemoglobin | HbO2 Soret + Q-bands, Hb, MetHb | 415, 542, 560, 577, 630 |
| Iron oxides (common) | Fe3+ hematite/goethite, Fe2+ pyroxene/olivine | 490, 530, 700, 900, 1000, 2000 |
| Iron oxides (expanded) | Ferrihydrite, lepidocrocite, siderite FeCO3 | 430, 475, 750, 1000, 1080, 2340 |
| Rare earths | Nd3+ (4 lines), Sm3+/Er3+, Pr3+ | 525, 580, 745, 800, 867, 940 |
| REE expanded | Dy3+, Yb3+, Er3+ (3 lines), Ho3+, Eu3+, Tb3+, Ce3+ | 400, 450, 490, 535, 537, 650, 910, 975, 1260, 1530 |
| Other metals | Cu2+, Cr3+, Co2+, Mn2+ | 435, 460, 550, 630, 700 |
| Water / ice | Liquid H2O, structural OH, ice | 970, 1200, 1380, 1450, 1940, 2000 |
| Clay minerals (common) | Al-OH, Mg-OH, NH4+ | 1410, 2160, 2200, 2250, 2310 |
| Clay minerals (expanded) | Pyrophyllite, gibbsite, boehmite, halloysite, phengite, palygorskite, sepiolite, talc | 1390, 2110, 2165, 2175, 2205, 2215, 2263, 2315, 2387, 2390 |
| Carbonates | CO3 calcite/dolomite/ankerite/siderite | 2120, 2320, 2330, 2335, 2350, 2500 |
| Sulphates (common) | Gypsum, alunite, jarosite | 1450, 1490, 1750, 2210, 2265 |
| Sulphates (expanded) | Epsomite, bassanite, kieserite, szomolnokite | 1000, 1490, 1640, 2165, 2200, 2220 |
| Organics | Cellulose, lignin, lipid, protein, starch | 1210, 1480, 1680, 1720, 2054, 2100 |
| Petroleum / coal | CH overtones, aromatic CH | 1210, 1720, 2300, 2310 |
| Polymers (common) | Polyethylene, polypropylene, nylon | 1150, 1215, 1660, 2310 |
| Polymers (expanded) | PET, PVC, polystyrene, rubber | 1640, 1680, 1720, 1730, 2030, 2170, 2240 |
| Silicates | Si-OH, Fe2+ feldspar | 1250, 1380, 2200 |
| Sheet silicates | Mica, amphibole, serpentine, chlorite, epidote | 2230, 2250, 2310, 2340, 2370 |
| Geological minerals | Prehnite, schorl/tourmaline, crocidolite | 440, 1490, 2120, 2243, 2320 |
| Urban materials | Portlandite Ca(OH)2 (hardened cement) | 1460 |
| Atmospheric refs | O2 A/B bands, H2O vapour | 690, 760, 820, 940, 1140 |
| Soil nutrients | Nitrate NO3-, ammonium NH4+ (2nd band) | 1410, 1560, 2050 |

## NOTES

### Input requirements

The module expects a 3D raster map created by *i.hyper.import* or any 3D
raster with wavelength metadata in one of the two supported formats:

- **r3.info history** — lines formatted as `Band N: WAVELENGTH nm, FWHM: F nm`
  (written by *i.hyper.import*)
- **r.support per-band metadata** — fields `wavelength=`, `FWHM=`, `valid=`,
  `unit=` in the description of `mapname#N` rasters

Spectra must be in surface reflectance units (0–1 range).  Apply atmospheric
correction (e.g., *i.hyper.atcorr* or *i.hyper.smac*) before using this module.

### Band depth computation

Band depth for band *i* at wavelength λ_i is:

```
depth_i = max(0, 1 − R_i / max(R_j  ∀j where |λ_j − λ_i| < 150 nm,
                                 j ≠ i, and j in the same coverage segment))
```

The continuum is the brightest neighbouring band within the ±150 nm window
**that the sensor actually measured next to the feature**. Most airborne and
drone products drop the two water-vapour windows, and a plain distance window
then reaches across the hole and takes the continuum from the far side, which
is not a shoulder of the feature at all. On an EMIT cube the gap at
1320-1439 nm is only 119 nm wide, so 10 of 244 bands were affected, and the
five bands carrying the 1400-1500 nm water feature read a depth of about 0.62
instead of 0.41.

A gap is detected relative to the sensor's own sampling: a spacing above
`max(50 nm, 5 × median spacing)`. That keeps a coarse multispectral sensor,
where 60 nm between bands is normal, in a single segment, while still
splitting a 5 nm imaging spectrometer at a 150 nm hole.

A band with no neighbour inside its own segment has no continuum and
therefore **no measurable depth**: it reports NaN in full-map mode and 0 in
point mode, not the `1 − R` that a default continuum of 1.0 used to produce.
That default invented a deep absorption at every isolated band, and the
deeper the darker the surface.

Typical thresholds:

| depth | Interpretation |
|-------|---------------|
| 0.00 | No absorption |
| 0.01–0.03 | Trace / below detection |
| 0.03–0.10 | Detectable — default threshold 0.02 |
| 0.10–0.30 | Moderate; most composite rules activate here |
| > 0.30 | Strong; nearly pure mineral or dense canopy |

### Storage scale of the reflectance values

Band depth is a ratio, so the feature detection does not care whether the cube
holds 0–1 reflectance or reflectance scaled by 10000, as surface-reflectance
products commonly do. The **absolute** thresholds do care:
**min_reflectance**, **max_reflectance** and **sand_albedo_threshold** are all
expressed in 0–1 reflectance.

**reflectance_scale=** is the divisor that converts stored values to
reflectance. With the default `0` it is read from the cube's own value range
and reported:

```
Reflectance scale detected: 10000 (stored value / 10000 = reflectance)
```

Detection picks 1, 100, 1000 or 10000 from the cube's maximum. Give the value
explicitly for a cube whose range does not reveal it — a scene that happens to
be entirely dark, for instance.

Note for existing command lines: before this, the absolute thresholds were
compared against the stored values directly, so working on a cube scaled by
10000 meant passing `sand_albedo_threshold=2500` to mean 0.25. That value now
means a reflectance of 2500 and will never be reached. Pass `0.25`.

### Detected, absent, or indeterminate

Band depth is a ratio, so it says nothing about the absolute size of the dip.
On a surface reflecting 2 % in the SWIR, a depth of 0.3 is a reflectance
difference of 0.006, at or below the noise of an atmospherically corrected
imaging spectrometer. Measured on the C2RD Pourrières ground targets, every
apparent polyamide detection — N-H 1535 nm and amide 2049 nm depths of
0.12 to 0.31 — came from a surface below 0.07 SWIR reflectance, and the
apparent depth tracked darkness rather than material.

A species whose own bands are all darker than **min_reflectance** is
therefore reported as **indeterminate**: neither detected nor absent. No rule
that needs that species, as a required or as a required-absent species, is
scored at such a pixel. Treating "unreadable" as "absent" is the same mistake
as treating an unsampled exclusion as a pass, and it has the same
consequence: the rule fires on exactly the surface that cannot support it.
With the default 0.02 this removes the open-water false positives on the
C2RD scene; raise it towards 0.07 for work on dark targets.

**validity=** writes the per-pixel outcome, which the class map alone cannot
express:

| code | meaning |
|---|---|
| NULL | the cube has no data here |
| 1 | classified |
| 2 | unclassified — readable, but no rule fired |
| 3 | low signal — every band below **min_reflectance** |
| 4 | saturated — a band above **max_reflectance** |

Low-signal and saturated pixels are left unclassified rather than passed off
as an ordinary no-match. On the C2RD scene the white glass-fibre panel reads
113 % reflectance at 550 nm with no data through the SWIR, and came out of
earlier runs merely "unclassified", indistinguishable from a confident
negative.

### Sensor band coverage

Before scoring, the module checks which of its composite rules this sensor can
actually decide, and reports the rest as warnings:

```
WARNING: Rule <kaolinite> is disabled for this sensor: required species
         <Al-OH_clay> is outside this sensor's bands.
```

A rule is disabled when the sensor samples none of the bands of one of its
required species, or none of the bands of one of its **required-absent**
species.  The second case matters more than it looks: an exclusion that cannot
be tested would otherwise pass everywhere, turning a veto into a no-op and
letting the rule fire on exactly the material it was written to reject.  Such
a rule is disabled rather than silently relaxed.

A band window counts as sampled only when a band falls near its nominal
centre, not merely somewhere inside the window.  Without that test a sensor
whose bands stop at 1482 nm would count as covering the 1400–1500 nm water
window, and the engine would read a band depth off the edge of a spectral
gap.  Cubes that drop the water-vapour windows — most airborne and drone
products do — are the normal case, not an exception.

### Per-band detection thresholds

A species' bands can differ in strength by a factor of two or more, so a
single species-level threshold serves none of them well.  A database entry may
therefore carry its own `min_depth`, which applies to that entry's bands
alone; the effective threshold for a band is the deepest of its own floor, the
rule's `min_depth_overrides` entry and the global **min_band_depth**.  The
open-water exclusion uses this: it vetoes on a depth of 0.60 at 970 nm, 0.55
at 1200 nm and 0.50 at 1940 nm, and is deliberately out of reach at 1450 nm,
where leaf water alone reaches 0.54 and so cannot separate water from
vegetation.

### Composite scoring

Each composite rule defines a `base_score` (0–1) and a list of required
species and optional bonus species.  A rule is scored at a pixel when every
required species is detected, no `required_absent` species is detected, and
every species the rule names is both sampled by this sensor and measurable
here. The reported confidence is then:

```
evidence_sp = sqrt(depth_margin_sp × band_support_sp)
evidence    = geometric mean over required species of evidence_sp
strength    = 0.55 + 0.45 × evidence
score       = min(1.0, base_score × strength × (1.0 + 0.05 × bonus_count))
```

- `depth_margin_sp` is how far the species' deepest band clears its own
  threshold, scaled over 0.20 of extra depth and clipped to [0, 1]
- `band_support_sp` is the fraction of that species' sampled bands which
  cleared the threshold
- `bonus_count` is the number of bonus species also detected

The `base_score` is therefore a ceiling that says how diagnostic the rule is
when it fires, and the pixel's own evidence decides how much of that ceiling
it earns. A match that only just clears every threshold keeps 55 % of the
base score; one well clear of all of them keeps all of it. Before this, the
confidence was `base_score × bonus_factor` and carried no information about
fit at all, which is how a strict run on the C2RD scene could report a mean
confidence of 0.99 over what was a degenerate two-way split.

The two components are deliberately independent: a single band scraping over
the threshold while its neighbours stay flat is a spike, whereas the same
depth across six of eight bands is an absorption band, and the geometric mean
requires both. Using the geometric mean across required species means one
species sitting at its threshold pulls the whole match down without zeroing
it, and adding a required species to a rule can never raise its score.

**Confidence is not a correctness probability.** It measures how well the
pixel fits the rule that was applied, which cannot detect a rule that is not
diagnostic in the first place. On the C2RD scene the surviving
misclassifications — a red roof, an orange court floor and a windsock all
read as `phycoerythrin_algae` — have genuinely deep features in the windows
the rule tests, and score high for that reason. Confidence calibration
narrows what a score means; it does not substitute for rule specificity.

### Band slice extraction

All Z-slices are extracted via `Rast3d_extract_z_slice()` from the GRASS
raster3d library (`libgrass_g3d`), which opens the 3D map with
`RASTER3D_NO_CACHE` and uses `Rast3d_get_block()` in tile-bulk mode.
Each tile at the target depth is loaded exactly once from disk, compared to
one function call per voxel in the default cache-based API.

### Atmospheric contamination flag

If bands matching O2 A/B-band (690, 760 nm) or H2O atmospheric (820, 940,
1140 nm) absorptions are detected above `min_band_depth`, a warning is issued
advising verification of atmospheric correction.

### SWIR overlap warning

When features are detected in the 2280–2390 nm zone and multiple competing
species (carbonate CO3, CH hydrocarbon, Mg-OH) are matched, a disambiguation
warning is issued.  Use FWHM as a discriminant: narrow features (< 30 nm)
indicate organic CH, broad features (> 80 nm) indicate minerals.

### Point mode output

In point mode (`-p`), the interpretation report is written to the console
via `gs.message()`.  The report lists:

- Dominant class and confidence
- All detected features above the depth threshold
- Chromophore assignments ranked by confidence
- Ranked material hypotheses (top 8 by default, all with `-v`)
- Any flag notes (atmospheric, overlap warnings)

### Flags

The **-n** flag restricts processing to only bands marked as `valid=1` in
metadata.  Useful for excluding atmospheric water vapour absorption bands
(around 1400 and 1900 nm) flagged during import.

The **-i** flag prints the band wavelength table and exits without processing
any raster data.  Use this first to verify sensor coverage and metadata
completeness.

The **-p** flag activates point mode.  A `coordinates=east,north` value is
required.  The spectrum at the specified location is read with a single
*r3.what* call and the full interpretation engine runs on that one pixel.
Extracting a 2D raster per band to read one pixel costs minutes on a large
cube; the per-band path is kept only as a fallback.

The **-v** flag prints the full hypothesis ranking in point mode.  Without
`-v`, only the top 8 hypotheses are shown.

### Output metadata

The class map (`output`) receives category labels automatically via
**r.category**, one label per class code.  Both output maps receive
descriptive metadata via **r.support** referencing the input map name.
The confidence map receives a `plasma` colour table via **r.colors**.

### Testing

The pure-Python interpretation engine can be tested without a running GRASS
session using the bundled `test_suite.py`:

```bash
python -m pytest test_suite.py -v
```

The suite covers **51 tests** across database integrity, `_band_depth()`,
`interpret_spectrum()` (edge cases and physics-based synthetic spectra for
green vegetation, liquid water, kaolinite, stressed vegetation, and new
classes), `_compute_depths_numpy()`, `_build_band_species_index()`,
`score_composites_numpy()`, chromophore confidence scaling, and targeted
tests for new classes (pyrophyllite, ferrihydrite, PVC, PET, gibbsite,
crocidolite, acid mine drainage, talc, phycoerythrin).

### Limitations

- Classification accuracy depends strongly on spectral coverage.  VNIR-only
  sensors miss clay, carbonate, and hydroxyl features (> 1400 nm).
- Spectra must be in surface reflectance (0–1).  Radiance or DN inputs will
  produce invalid results.
- The ±150 nm continuum window may blend features that are closely spaced
  (e.g., the 2330/2350 nm carbonate doublet is partially merged).
- Dense vegetation mixtures with strong water absorption may trigger multiple
  competing hypotheses with similar scores.
- The chrysotile\_asbestos class (27) requires both chrysotile Mg-OH (2320 nm)
  and antigorite Mg-OH (2370 nm) features; low sensor SWIR2 SNR may miss it.

## EXAMPLES

::: code

    # Check band wavelengths and metadata before processing
    i.hyper.spectroscopy input=enmap_boa -i

:::

::: code

    # Full-map classification: dominant class + confidence
    i.hyper.spectroscopy input=enmap_boa \
                         output=enmap_class \
                         confidence=enmap_conf

:::

::: code

    # Use only valid bands, lower detection threshold for subtle features
    i.hyper.spectroscopy input=enmap_boa \
                         output=enmap_class \
                         confidence=enmap_conf \
                         min_band_depth=0.01 \
                         -n

:::

::: code

    # Point-mode interpretation at a known mineral outcrop
    i.hyper.spectroscopy input=enmap_boa \
                         coordinates=450320,5317880 \
                         -p -v

    # Example console output:
    # ──────────────────────────────────────────────────────────────────────
    #   SPECTRAL INTERPRETATION REPORT
    # ──────────────────────────────────────────────────────────────────────
    #   Dominant class : kaolinite
    #   Active features: 7
    #   Chromophores   : 5
    #   Hypotheses     : 4
    #
    #   RANKED HYPOTHESES (top 8)
    #   [ 1] kaolinite                           mineral      0.88  █████████████████
    #   [ 2] muscovite_illite                    mineral      0.78  ███████████████
    #   [ 3] advanced_argillic_alteration        geological   0.74  ██████████████
    #   [ 4] atomic_Al-OH                        clay         0.49  █████████
    #
    #   CHROMOPHORE DETECTIONS
    #   Al-OH_clay                          conf=0.72  [1412nm]
    #   Al-OH                               conf=0.68  [2203nm]
    #   structural_OH                       conf=0.44  [1389nm]
    # ──────────────────────────────────────────────────────────────────────

:::

::: code

    # Typical full workflow from import to classification
    # 1. Import
    i.hyper.import input=/data/ENMAP_L2A.hdf product=enmap output=enmap

    # 2. Atmospheric correction
    i.hyper.smac input=enmap output=enmap_boa

    # 3. Continuum removal (optional — sharpens absorption depths)
    i.hyper.continuum input=enmap_boa output=enmap_cr

    # 4. Classify
    i.hyper.spectroscopy input=enmap_cr \
                         output=enmap_class \
                         confidence=enmap_conf \
                         -n

    # 5. Visualise
    d.rast map=enmap_class
    d.legend map=enmap_class
    r.univar map=enmap_conf zones=enmap_class

:::

::: code

    # Extract all pixels classified as chrysotile asbestos (class 27)
    r.mapcalc expression="asbestos_mask = if(enmap_class == 27, 1, null())"
    r.colors map=asbestos_mask color=reds
    r.stats -a input=asbestos_mask

:::

::: code

    # Combine spectroscopy with geology for alteration mapping
    i.hyper.geology input=enmap_cr \
                    output_alteration=enmap_alteration

    # Pixels where spectroscopy sees kaolinite (6) AND geology sees argillic (4)
    r.mapcalc expression="argillic_confirmed = \
        if(enmap_class == 6 && enmap_alteration == 4, 1, null())"

:::

## SEE ALSO

[i.hyper.atcorr](i.hyper.atcorr.html),
[i.hyper.continuum](i.hyper.continuum.html),
[i.hyper.geology](i.hyper.geology.html),
[i.hyper.import](i.hyper.import.html),
[i.hyper.indices](i.hyper.indices.html),
[i.hyper.smac](i.hyper.smac.html),
[r.category](https://grass.osgeo.org/grass-stable/manuals/r.category.html),
[r.colors](https://grass.osgeo.org/grass-stable/manuals/r.colors.html),
[r.mapcalc](https://grass.osgeo.org/grass-stable/manuals/r.mapcalc.html),
[r.univar](https://grass.osgeo.org/grass-stable/manuals/r.univar.html),
[r3.info](https://grass.osgeo.org/grass-stable/manuals/r3.info.html),
[r3.what](https://grass.osgeo.org/grass-stable/manuals/r3.what.html)

## REFERENCES

- Hunt, G. R. (1977). Spectral signatures of particulate minerals in the
  visible and near infrared. *Geophysics*, 42(3), 501–513.
- Clark, R. N., et al. (1990). High spectral resolution reflectance
  spectroscopy of minerals. *Journal of Geophysical Research: Solid Earth*,
  95(B8), 12653–12680.
- Kokaly, R. F., et al. (2017). USGS Spectral Library Version 7.
  U.S. Geological Survey Data Series 1035.
  https://doi.org/10.3133/ds1035
- Clevers, J. G. P. W., & Kooistra, L. (2012). Using hyperspectral remote
  sensing data for retrieving canopy chlorophyll and nitrogen content.
  *IEEE Journal of Selected Topics in Applied Earth Observations and Remote
  Sensing*, 5(2), 574–583.
- Vane, G., & Goetz, A. F. H. (1993). Terrestrial imaging spectrometry:
  current status, future trends. *Remote Sensing of Environment*, 44(2–3),
  117–126.
- Koerting, F., et al. (2024). VNIR-SWIR Imaging Spectroscopy for Mining:
  Insights for Hyperspectral Drone Applications. *Mining*, 4, 1013–1057.
  https://doi.org/10.3390/mining4040057

## AUTHORS

Created for the i.hyper module family (2025).
