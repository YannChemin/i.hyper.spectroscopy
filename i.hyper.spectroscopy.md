## DESCRIPTION

*i.hyper.spectroscopy* performs per-pixel spectral feature interpretation
of hyperspectral imagery imported as 3D raster maps (`raster_3d`) by
[i.hyper.import](i.hyper.import.html).

The module reads wavelength metadata from hyperspectral 3D raster bands,
computes local band depths across the VNIR-SWIR range (400–2500 nm), matches
detected absorption features to a physics-based database of 70+ absorbing
species, and scores 28 composite material hypotheses per pixel.  Two output
rasters are produced: a dominant-class integer map and a confidence map.

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

### Absorption feature database

The built-in database covers 70+ entries across the full VNIR-SWIR range:

| Group | Key species | Diagnostic wavelengths (nm) |
|-------|-------------|------------------------------|
| Vegetation pigments | Chlorophyll a+b, carotenoids, phycocyanin | 430, 480, 620, 665, 710 |
| Hemoglobin | HbO2 Soret + Q-bands, Hb, MetHb | 415, 542, 560, 577, 630 |
| Iron oxides | Fe3+ hematite/goethite, Fe2+ pyroxene/olivine | 490, 530, 700, 900, 1000, 2000 |
| Rare earths | Nd3+ (4 lines), Sm3+/Er3+, Pr3+ | 525, 580, 745, 800, 867, 940 |
| Other metals | Cu2+, Cr3+, Co2+, Mn2+ | 435, 460, 550, 630, 700 |
| Water / ice | Liquid H2O, structural OH, ice | 970, 1200, 1380, 1450, 1940, 2000 |
| Clay minerals | Al-OH, Mg-OH, NH4+ | 1410, 2160, 2200, 2250, 2310 |
| Carbonates | CO3 (calcite, dolomite) | 2120, 2330, 2350, 2500 |
| Sulphates | Gypsum, alunite, jarosite | 1450, 1490, 1750, 2210, 2265 |
| Organics | Cellulose, lignin, lipid, protein, starch | 1210, 1480, 1680, 1720, 2054, 2100 |
| Petroleum / coal | CH overtones, aromatic CH | 1210, 1720, 2300, 2310 |
| Polymers | Polyethylene, polypropylene, nylon | 1150, 1215, 1660, 2310 |
| Silicates | Si-OH, Fe2+ feldspar | 1250, 1380, 2200 |
| Sheet silicates | Mica, amphibole, serpentine, chlorite, epidote | 2230, 2250, 2310, 2340, 2370 |
| Atmospheric refs | O2 A/B bands, H2O vapour | 690, 760, 820, 940, 1140 |

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
depth_i = max(0, 1 − R_i / max(R_j  ∀j where |λ_j − λ_i| < 150 nm, j ≠ i))
```

The continuum is the brightest neighbouring band within the ±150 nm window.
When no neighbour exists (isolated band), continuum = 1.0.  Bands with
`depth < min_band_depth` are not counted as absorbed.

Typical thresholds:

| depth | Interpretation |
|-------|---------------|
| 0.00 | No absorption |
| 0.01–0.03 | Trace / below detection |
| 0.03–0.10 | Detectable — default threshold 0.02 |
| 0.10–0.30 | Moderate; most composite rules activate here |
| > 0.30 | Strong; nearly pure mineral or dense canopy |

### Composite scoring

Each composite rule defines a `base_score` (0–1) and a list of required
species and optional bonus species.  When all required species are detected:

```
score = min(1.0, base_score × (1.0 + 0.05 × bonus_count))
```

where `bonus_count` is the number of bonus species that are also detected.
Rules with `required_absent` species fail if any of those species are
detected at the pixel.

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
required.  The spectrum at the specified location is extracted (one call to
`r.what` per band) and the full interpretation engine runs on that single pixel.

The **-v** flag prints the full hypothesis ranking in point mode.  Without
`-v`, only the top 8 hypotheses are shown.

### Output metadata

The class map (`output`) receives category labels automatically via
**r.category**, one label per class code.  Both output maps receive
descriptive metadata via **r.support** referencing the input map name.
The confidence map receives a `plasma` colour table via **r.colors**.

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

[i.hyper.import](i.hyper.import.html),
[i.hyper.geology](i.hyper.geology.html),
[i.hyper.continuum](i.hyper.continuum.html),
[i.hyper.indices](i.hyper.indices.html),
[i.hyper.atcorr](i.hyper.atcorr.html),
[i.hyper.smac](i.hyper.smac.html),
[r.mapcalc](https://grass.osgeo.org/grass-stable/manuals/r.mapcalc.html),
[r.category](https://grass.osgeo.org/grass-stable/manuals/r.category.html),
[r.colors](https://grass.osgeo.org/grass-stable/manuals/r.colors.html),
[r.univar](https://grass.osgeo.org/grass-stable/manuals/r.univar.html),
[r3.info](https://grass.osgeo.org/grass-stable/manuals/r3.info.html)

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
