# i.hyper.spectroscopy — planned improvements

Priorities below come from the C2RD Pourrieres drone campaign (HySpex
PARGE/ATCOR mosaic, 8471 x 7284 x 408 at 3 cm, 411-2450 nm), where the module
was run against 23 surveyed ground targets and 18 surveyed infrastructure
classes. Each item states the evidence that motivated it and what would count
as done, so none of them has to be re-derived.

The scene matters because it is the normal case, not a corner case: like most
airborne and drone products it **drops both water-vapour windows**, so it has
no bands in 1329.0-1482.6 nm nor in 1784.5-1984.0 nm. Several items below
exist only because the engine assumed contiguous coverage.

## Done

- **Per-band detection thresholds.** A feature-database entry may carry its own
  `min_depth`, applying to that entry's bands alone. The effective threshold
  for a band is the deepest of its own floor, the rule's `min_depth_overrides`
  entry and the global `min_band_depth`. A species whose bands differ in
  strength by more than a factor of two (the water OH overtones) could not be
  expressed before.
- **Open-water exclusion rebuilt as `water_body_OH_exclusion`.** Four
  overtones with measured floors: 0.60 at 970 nm and 0.55 at 1200 nm (on this
  scene open water reaches 0.80-0.89 and 0.75-0.78, the deepest non-water
  surface 0.35 and 0.31), 0.50 at 1940 nm carried over from the laboratory
  branch and unvalidated, and 1450 nm deliberately out of reach at 0.80
  because leaf water alone reaches 0.54 there and would make it a vegetation
  veto rather than a water veto. It replaces `water_body_OH_1940`, which this
  sensor does not sample at all.
- **Sensor coverage check.** `unsupported_rules()` disables any rule whose
  required species, or whose required-absent species, has no band near its
  nominal centre in this sensor, and says so. Previously an untestable
  exclusion passed everywhere, turning a veto into a no-op: that is why the
  polymer-sheet rule fired on 587 m2 of pond and pool water against roughly
  4 m2 of real sheeting. Six rules on this cube were failing silently.
- **All five camouflage rules now carry the water veto.** Only the field sheet
  rule had one; the three laboratory fabric rules keyed on C-H overtones
  alone, and deep open water satisfied them.
- **Point mode reads the spectrum with `r3.what`.** See item 1.

## 1. Point mode via r3.what — done, keep as the pattern

Point mode used to call `extract_band()` plus `r.what` once per band: it
materialised one full-scene 2D raster per band to read a single pixel. On this
cube that is 408 full-scene rasters for one spectrum. `r3.what` reads only the
voxels asked for and returns the whole column in one call. The slicing loop is
kept as a fallback for when `r3.what` is unavailable or returns an unexpected
band count.

The same lever applies elsewhere: for point and polygon sampling, `r3.what`
beat `r3.to.rast` slicing by about four orders of magnitude in the campaign
work (16,482 points x 408 bands in 2.1 s). Any future per-point or per-polygon
extraction path should use it.

### The same change, applied to the three sibling modules — done

An audit of every `i.hyper.*` module found exactly four with a point mode, and
three carried a verbatim copy of the loop replaced here. All three now take
`_spectrum_via_r3what()` with the per-band slicing loop kept as a fallback:

| Module | Function | Before (measured) | After |
|---|---|---|---|
| *i.hyper.cannabis* | `read_pixel_spectrum()` | 5.54 s/band, ~37.7 min for 408 | 0.60 s |
| *i.hyper.poppy* | `read_pixel_spectrum()` | as above | 0.61 s |
| *i.hyper.sleuth* | `read_pixel_spectrum()` | as above | 0.60 s |

All four modules return identical values at the same coordinate, and the
helper indexes by band number rather than by position, so `only_valid`
filtering and any wavelength sorting stay correct.

One pre-existing defect surfaced while building: *i.hyper.poppy* ships only
`i.hyper.poppy.html` and no `i.hyper.poppy.md`, so `make` fails at the man
page step (`No rule to make target .../i.hyper.poppy.md`). The script itself
installs; only the manual does not build.

Every other `extract_band()` / `r3.to.rast` caller in the family is **not** a
candidate: *i.hyper.albedo*, *i.hyper.cropname*, *i.hyper.endmembers*,
*i.hyper.geology*, *i.hyper.methane*, *i.hyper.rgb*, *i.hyper.smac* (and its
`lib/aod.py`, `lib/o3.py`, `lib/utils.py`, `lib/wvc.py`) and
*i.hyper.traficability* all slice because they need every pixel of a band as a
2D map. `r3.what` would be the wrong tool there: it reads only the voxels
asked for, which is the whole point for a point query and a pessimisation for
a full raster. *i.hyper.smac* already has the better answer for that case in
`lib/r3d_slice.py`, which reads a Z-slice through ctypes with no subprocess at
all.

## 2. Continua must respect spectral gaps — correctness bug

`_compute_depths_numpy()` and `_band_depth()` take the continuum as the
brightest band within +-150 nm of the feature, with no regard for whether the
sensor actually covers that interval. Across the 1784-1984 nm gap the nearest
"neighbour" is 200 nm away, so every band depth near a gap edge is measured
against the far side of the gap. This is how a water band depth of 0.64-0.70
came to be reported for a cube that has no bands within 150 nm of 1940 nm.

- Split the band list wherever the spacing exceeds about 50 nm, and compute
  continua within each contiguous segment only.
- Report no depth, rather than a wrong one, where a feature has no continuum
  shoulder on at least one side inside its own segment.
- Acceptance: on a cube with a 200 nm gap, no feature within 150 nm of a gap
  edge takes its continuum from the other side.

## 3. Three-valued species logic: detected / absent / indeterminate

The water bug was `absent` and `untestable` being the same value. That is now
fixed at rule granularity, but the same confusion remains per pixel.

All six apparent polyamide detections among the field targets came from
surfaces with SWIR reflectance below 700/10000, where band depth is noise: the
apparent amide depths of 0.12-0.31 correlate monotonically with darkness, not
with material. At the same time the eight bright sheets (24-35 % SWIR
reflectance) show genuinely no C-H absorption, so absence there is real signal.
The engine cannot currently tell those two situations apart.

- Add a minimum-reflectance guard per feature, below which a species is
  *indeterminate* rather than detected or absent.
- A rule with an indeterminate input declines to fire rather than firing
  positive.
- Saturated pixels are a second case of the same kind: the white glass-fibre
  panel reads 113 % reflectance at 550 nm with nodata through SWIR1 and comes
  out merely "unclassified", indistinguishable from a confident negative.
- Acceptance: a validity or quality output map alongside class and confidence,
  distinguishing classified, unclassified, low-signal and saturated.

## 4. A shape detector beside the feature detector — architectural

The engine can express only discrete absorption minima. The campaign's polymer
sheets have none: for one white sheet it finds 8 features in 408 bands and
nothing at all between 2050 and 2300 nm. Their signature is broad continuum
curvature, so **no rule in the present formalism can state it**, at any
threshold. A field-derived 2122/2245 nm pair reached 11/12 recall when measured
with a sliding-window continuum, and still fired on nothing once wired in,
because the engine's test is a different measurement.

- Add band-ratio and continuum-slope descriptors between anchor wavelengths as
  a second species kind, so a rule can require a shape and not only a minimum.
- Acceptance: a rule that detects the 12 surveyed sheets without firing on
  bare soil, limestone or vegetation in the same scene.

## 5. A `library=` option for spectral-angle or matched-filter matching

The highest-value single addition. Region growing at SAM < 0.20 rad already
delineated all 12 surveyed sheets cleanly, mean within-region angle
0.047-0.148, where the rule engine detected none of them. The library file
format already exists from the campaign: wavelength column plus per-class mean
and standard deviation.

- Read a labelled library CSV, resample it to the cube's bands, and output
  best-match class plus the angle to it.
- Treat the angle as the confidence, which unlike the rule score has a
  physical meaning.
- Acceptance: on the C2RD mosaic, the 12 sheets are recovered from
  `c2rd_target_spectra.csv` with an angle below 0.20 rad, and open water and
  limestone are not.

## 6. Calibrated confidence

`score = base_score x (1 + 0.05 x bonus_count)` carries no information about
how strong the evidence actually is. A strict run on this scene reached mean
confidence 0.99, minimum 0.90, on what was a degenerate five-class split
driven entirely by the NDRE pre-screen.

- Fold depth margin above threshold, and the number of supporting bands, into
  the reported confidence.
- Acceptance: confidence correlates with agreement against the surveyed ground
  classes, instead of being near-constant.

## 7. Regression testing against field data, not synthetic spectra

`test_suite.py` builds Gaussian bumps. Every false positive found in the
campaign was invisible to it, and the surveyed data is itself a confusion
matrix: three independent red and orange man-made surfaces (a roof, a court
floor, a windsock) all classify as `phycoerythrin_algae`; white paint and
court strips fall through to the quartzofeldspathic sand fallback; deep pool
water classifies as `acid_mine_drainage`; the genuinely eutrophic pond
classifies as `ice_snow`; and asphalt is missed entirely despite
`asphalt_bitumen` existing as a class.

- Add `tools/validate_against_library.py` taking a labelled library CSV and
  printing per-class precision and recall.
- Wire the 23 surveyed targets and 18 infrastructure classes in as a gate, so
  a knowledge-base edit is tested against real ground.
- Acceptance: editing a rule that changes behaviour on surveyed ground fails
  the gate rather than passing unnoticed.

## 8. Smaller items

- A print mode that lists, for a given cube, which rules this sensor can
  decide and which it cannot, so coverage is inspectable before a long run
  rather than reported as warnings during one.
- Generic plastic rules (`PVC_plastic`, `PET_plastic`,
  `polyethylene_plastic`, `polystyrene`, `rubber_polyisoprene`) fire on
  almost everything in point mode at the default `min_band_depth=0.02`,
  including asphalt and bare soil. They need the same treatment as the
  camouflage branch: floors measured against real surfaces.
- Class 0 is exported as NULL, so "unclassified" and "no data" cannot be
  distinguished in the output map.
- The module leaves the 3D computational region changed on exit; it should
  restore it through a scoped region manager.
- Raster writing round-trips through ASCII at `%.4f`, which both costs time
  and silently quantises reflectance.
