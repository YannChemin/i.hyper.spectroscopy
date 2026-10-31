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
- **Continua respect spectral gaps**, and a band with no shoulder reports no
  depth instead of `1 - R`. See item 2.
- **Three-valued species logic** (detected / absent / indeterminate), plus a
  `validity=` output map. See item 3.
- **Calibrated confidence** from depth margin and band support. See item 6.
- **Reflectance storage scale is resolved explicitly** (`reflectance_scale=`,
  detected from the cube's range by default). Found while testing item 3: band
  depth is a ratio and never cared, but every absolute threshold did, and on a
  cube scaled by 10000 the new saturation guard flagged all 346,800 pixels of
  a test window at once. This is also why `sand_albedo_threshold` previously
  had to be given as `2500` to mean 0.25; it now means what it says, which is
  a **breaking change for existing command lines**.

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

Two pre-existing defects surfaced alongside and are now fixed: *i.hyper.poppy*
shipped only `i.hyper.poppy.html` and no `i.hyper.poppy.md`, so `make` failed
at the man page step; and *i.hyper.cannabis* and *i.hyper.poppy* had fallen
behind *i.hyper.sleuth* on two shared-code improvements (`garray.array(mapname=)`
instead of the two-step read, and in-place `_safe_cube()` via `np.nan_to_num`).
All three are in sync again.

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

## 2. Continua must respect spectral gaps — done

`_spectral_segments()` numbers the contiguous blocks of coverage, splitting
wherever the band spacing exceeds `max(50 nm, 5 x median spacing)`, and both
`_compute_depths_numpy()` and `_band_depth()` now take the continuum only from
bands in the feature's own segment. The relative threshold is what makes this
safe for coarse sensors: 60 nm between bands is normal for a multispectral
instrument and must not read as a gap.

A band with no neighbour inside its own segment has no continuum and so no
measurable depth. It now reports NaN (full map) or 0 (point mode) instead of
the `1 - R` that a default continuum of 1.0 produced, which invented a deep
absorption at every isolated band, deeper the darker the surface.

Measured on a real EMIT cube (`emit_20250410`, 244 bands): the gap at
1320-1439 nm is only 119 nm wide, narrower than the +-150 nm window, so **10
of 244 bands** were taking their continuum from the far side. On 64 sampled
pixels the five bands carrying the 1400-1500 nm water feature read a mean
depth of 0.62 where the correct value is 0.41, a mean error of +0.21 and a
maximum of +0.34. Every unaffected band is bit-identical before and after.

Three of the module's own tests (`ferrihydrite`, `crocidolite_not_chrysotile`,
`acid_mine_drainage`) were passing *because* of this bug: their synthetic
spectra place the nearest continuum band 170 nm from the feature, so the
features had no shoulder and the 1.0 fallback manufactured the depth. Their
spectra now carry a continuum anchor within 150 nm of each diagnostic band.

## 3. Three-valued species logic — done

`_species_measurable()` marks a species *indeterminate* at a pixel when all of
its own bands are darker than **min_reflectance** (default 0.02 reflectance).
No rule is scored where any of its required or its required-absent species is
indeterminate, in both the full-map and the point path: "unreadable" is not
"absent", and treating it as such is the same mistake as treating an unsampled
exclusion as a pass.

`validity=` writes the per-pixel outcome the class map cannot express: NULL
for no data, 1 classified, 2 unclassified (readable, no rule fired), 3 low
signal, 4 saturated. Low-signal and saturated pixels are left unclassified
rather than passed off as an ordinary no-match.

Verified on the surveyed spectra: at the default 0.02 the guard marks 85 to 90
species indeterminate on pool and pond water and **removes the
`acid_mine_drainage` false positive on open water** that the POI cross-tab
found. A unit test holds the same spectrum at full brightness and at a
hundredth of it: the band depths are identical to 1e-9, the bright pixel
fires, the dark one is withheld, and with the guard disabled the dark one
fires again.

What it does not reach: the dark-target amide artefacts at 0.023-0.034
reflectance need a floor near 0.07 to be excluded, and 0.07 is too aggressive
as a default because a vegetation pixel's chlorophyll bands sit at 0.03-0.08
and would be blocked with it. The option exists for dark-target work; the
default is set where it demonstrably removes false positives without
suppressing real classes.

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

## 5. A `library=` option, but in `i.hyper.sleuth`, not here

Still the highest-value single addition, and smaller than first written.
Region growing at SAM < 0.20 rad already delineated all 13 surveyed sheets
cleanly, mean within-polygon angle 0.046-0.105, where the rule engine detected
none of them. The library file format already exists from the campaign:
wavelength column plus per-class mean and standard deviation.

**The matching mathematics is already implemented.** `i.hyper.sleuth` offers
nineteen measures -- `sam`, `sid`, `sid_sam`, `ed`, `sad`, `sca`, `cr_sam`,
`cr_ed`, `gd1`, `gd2`, `xcorr`, `dtw`, `ssim`, `jsd`, `bhatt`, `mtf`, `cem`,
plus `ensemble` and `consensus` with six fusion modes, rank calibration and
diversity weighting. Writing a second SAM in this module would duplicate it
and be worse. What sleuth lacks is **multi-class**: `reference=` and
`reference_file=` each take one spectrum and `output=` is one similarity
raster.

So the work belongs in `i.hyper.sleuth`:

- Add `library=`, reading a labelled CSV whose columns are
  `wavelength_nm`, then `<class>_mean` and `<class>_sd` per class -- the shape
  `c2rd_target_spectra.csv` already has.
- Resample each class to the cube's bands through the existing
  `WavelengthLUT`, loop the existing matcher over classes, and output
  best-match class plus the angle to it (`output=` class map,
  `confidence=` angle).
- Treat the angle as the confidence, which unlike the rule score has a
  physical meaning. The existing `consensus` machinery then comes for free.
- Keep `min_reflectance`-style guards out of it: an angle is scale-free, which
  is most of the reason it works where the depth rules do not.
- Acceptance: on the C2RD mosaic, the 13 sheets are recovered from
  `c2rd_target_spectra.csv` with an angle below 0.20 rad, and open water and
  limestone are not.
- Then this module's `library=` is just a delegation, or nothing at all.

## 6. Calibrated confidence — done, with a limit worth stating

The reported confidence is now

```
evidence_sp = sqrt(depth_margin_sp x band_support_sp)
evidence    = geometric mean over required species of evidence_sp
score       = min(1, base_score x (0.55 + 0.45 x evidence) x (1 + 0.05 x bonus))
```

where `depth_margin_sp` is how far the species' deepest band clears its own
threshold over a 0.20 scale, and `band_support_sp` is the fraction of that
species' sampled bands which cleared it. The base score became a ceiling that
says how diagnostic the rule is; the pixel's own evidence decides how much of
it is earned. A match that only just clears every threshold keeps 55 % of the
base score, which is what a unit test checks.

The two components are independent on purpose: one band scraping over the
threshold while its neighbours stay flat is a spike, the same depth across six
of eight bands is an absorption band, and the geometric mean needs both. The
geometric mean across species means one species at its threshold pulls the
match down without zeroing it, and adding a required species can never raise a
rule's score.

**The acceptance criterion as written was not met, and cannot be met this
way.** Confidence does not correlate with agreement against the surveyed
classes, because the surviving errors are not weak-evidence errors: the red
roof, the orange court floor and the windsock that all read
`phycoerythrin_algae` have genuinely deep features in the windows that rule
tests, and earn high confidence honestly. Across 15 surveyed spectra the top
hypothesis now takes 8 distinct confidence values instead of 3, spanning
0.78-0.99, but the wrong answers are not the low ones. Calibration narrows
what a score means; it cannot detect a rule that is not diagnostic. That is
item 4 and item 5.

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
- Class 0 is still exported as NULL in the class map, so "unclassified" and
  "no data" remain indistinguishable *there*; the `validity=` map now carries
  the distinction, but the class map itself should stop conflating them.
- The full-map path leaves the 3D computational region changed on exit; it
  should restore it through a scoped region manager, as the point-mode
  `r3.what` helper now does with `gs.use_temp_region()`. A run killed
  mid-flight leaves the region at a single depth, which then silently
  changes what a later `r3.what` returns.
- Raster writing round-trips through ASCII at `%.4f`, which both costs time
  and silently quantises reflectance.
- The `atexit` cleanup of temporary 2D rasters does not survive SIGTERM: a
  run stopped part-way leaves one full-scene raster per band behind (143
  rasters, 254 MB, after one interrupted point-mode run on the C2RD cube).
  A signal handler should remove them.
- The sand fallback's `sand_albedo_threshold` and the two reflectance guards
  are the only scale-dependent numbers in the module. Any further absolute
  threshold belongs behind `reflectance_scale` too; every other test in the
  engine is a ratio and is indifferent to storage scale.
- `min_reflectance` guards a species by the brightest of its own bands. A
  contrast-based test would be stricter and better founded -- a dip matters
  when `continuum x depth` exceeds the sensor noise, not when the surface is
  bright -- but it needs the continuum value, which the depth array does not
  carry. Worth doing when the depth computation is next touched: it would
  catch the dark-target amide artefacts at 0.023-0.034 reflectance that the
  present default cannot reach.
- `reflectance_scale` detection reads the cube's maximum, so a scene that is
  entirely dark could be misread. It is reported on every run and can be
  given explicitly, but a check against the band metadata's own scale factor
  (where `i.hyper.import` records one) would be better than inference.

## 9. Cross-module items found while testing the siblings

These came out of writing `i.hyper.poppy`'s unit suite and running the
`test_ras3d_*` suites against the real Wyvern and Tanager cubes. They are
recorded here because this is the family's working TODO, but four of them are
not about `i.hyper.spectroscopy`.

**Fixed.**

- `_convert_wl_nm` accepted `micrometer` and `micron` but not
  **`micrometers`**, which is the spelling ENVI headers use and which
  `i.hyper.import/i_hyper_lib/envi.py` already accepts. A cube imported from
  an ENVI header had 1.45 um read as 1.45 nm, with only a warning -- a silent
  1000x error. The docstrings in `sleuth` and `cannabis` even claimed
  `micrometer(s)` was accepted. Fixed in `poppy`, `sleuth`, `cannabis` and
  `spectroscopy`, now taking the plural and British spellings too, with 20
  parametrized cases pinning it.
- `GLOBAL_STATS_METHODS = {'mtf', 'cem'}` was defined in `poppy`, `sleuth`
  and `cannabis` and consulted by nothing in point mode, so both methods ran
  on a single pixel, where the scene covariance they need is undefined
  (`n - 1 = 0`). `match_mtf` built an all-NaN covariance and only degraded
  gracefully by accident, because `pinv` happened to raise and hit the SAM
  fallback; `match_cem` returned a degenerate score silently. Point mode now
  skips them with a warning in all three.
- The `test_ras3d_*` suites called `importlib.import_module('i_hyper_<x>')`
  while the file on disk is `i.hyper.<x>.py`, which is not an importable name
  however `sys.path` is set. **Nine suites had therefore never once imported
  their module.** `cropname` and `rgb` already did it correctly with
  `spec_from_file_location`; that is now a shared `load_addon()` helper in
  `test_ras3d_common.py`, distributed to all 12 copies.

**Open.**

- `match_sam` adds `EPS` inside both norms instead of guarding the division,
  so a perfect match returns `cos(a) = 1 - O(EPS)` and `arccos` amplifies that
  to `sqrt(EPS)`, about 1e-6. Bounded and cosmetic, at the one point where the
  answer is already unambiguous. Left alone rather than touch four modules'
  matcher core; the poppy suite documents it with a 1e-5 tolerance and says
  why.
- **Three `.wl.json` readers are dead code, and three modules have none.**
  The canonical writer is `ras3d_write` (`/usr/lib/python3/dist-packages/ras3d_write.py`,
  lines 104-107): the path is the cube path with its **extension stripped**
  plus `.wl.json`, and the content is a **flat JSON list of floats**,
  `json.dump([round(float(w), 6) for w in wavelengths], f)`. Measured against
  that contract:
  - **Correct**: `albedo`, `cropname`, `geology`, `poppy`, `rgb` -- they try
    the suffix list, read a flat list and scale um to nm.
  - **Broken, and could never have worked**: `cannabis`, `sleuth`,
    `spectroscopy`. All three build the path as `raster3d + '.wl.json'`, which
    keeps the `.tiff` and so never finds the file, *and* iterate it as
    `[{'band': .., 'wavelength': .., 'fwhm': .., 'valid': ..}, ..]`, a schema
    nothing writes. Doubly dead: the wrong path means the wrong schema was
    never even reached. Fix all three to the canonical contract, keeping their
    own `min_wl`/`max_wl`/`only_valid` filtering, with `albedo`'s block as the
    model.
  - **Absent**: `continuum`, `indices`, `specresamp` resolve wavelengths only
    through `$MAPSET/grid3/<name>/hyper.json` and fatal on a cube opened
    straight from a file. Adding the fallback is a small feature, six lines
    copied from `albedo`, but it is a feature rather than a test fix.

  The `test_wavelength_sidecar` test in each `test_ras3d_*` suite was written
  to check exactly this and had never run, because of the import bug above.
