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

## Implementation order

Ranked by architectural benefit -- how much each one changes what the rest of
the family can rely on -- rather than by raw value or by effort. Section
numbers in brackets point at the full entry.

**Shared contracts that are currently wrong.** These are inherited by every
caller, so they come first whatever else is pending.

1. **Gaps as no-support in the lookup service** [10, first bullet] --
   **done**, and it moved 11 of 20 C2RD identifications. It also showed why
   a relative overlap threshold alone is not enough: the denominator is the
   query's own coverage, so the floor loosens as the query gets gappier. The
   service now takes an absolute `min_overlap_bands` too. The rest of the
   family still has to pass one: every module that calls the service with
   `min_overlap_fraction=` alone is measuring against a moving denominator.
2. **The `.wl.json` canonical contract across six modules** [9]. One writer,
   six divergent readers: three resolve a path that cannot exist *and*
   iterate a schema nothing writes, three have no reader and fatal on a cube
   opened from a file. A contract exists and no implementation matches it.
3. **A shared service client in `i_hyper_lib`** [10]. Do it before more
   callers arrive: `i.hyper.sleuth` becomes the third consumer at item 6.
   One helper also gets graceful degradation when the service is down,
   instead of each module failing its own way.

**Makes everything after it verifiable.**

4. **Field-data regression gate** [7]. Until the surveyed ground is a gate,
   items 6, 7 and 9 cannot be shown to have helped -- which is why this
   outranks them despite being plumbing.
5. **Region and temp-raster lifecycle** [8, two bullets]. An interrupted run
   leaves the 3D region changed, silently altering what a later `r3.what`
   returns, and `atexit` does not survive SIGTERM. Small work, removes a
   whole class of silent failure.

**New capability, once the ground is solid.**

6. **Multi-class `library=` in `i.hyper.sleuth`** [5]. Consolidates matching
   where the nineteen measures already live, so no second SAM is written
   here and this module's `library=` becomes a delegation or nothing. Also
   step 2 of the extraction pipeline.
7. **Shape detector as a second species kind** [4]. Broad continuum
   curvature cannot be stated in the present formalism at any threshold.
   Needs item 4 first or its acceptance test is unfalsifiable.
8. **Drop the ASCII round-trip for raster writing** [8]. `%.4f` through
   `r.in.ascii` costs time and silently quantises reflectance; the I/O path
   is wrong for the data type.

**High value, low architecture.**

9. **Category filtering as documented practice** [9]. The cheapest real
   improvement on this list: the options exist and the server applies them,
   so the work is deciding which categories a field scene is restricted to
   and measuring before/after. Pair it with the `homogeneity=` run [9] as
   its validation.
10. **Service packaging as a systemd user unit** [10]. Worth little before
    item 3 and quite a lot after it: it turns the service from something
    started by hand into a dependency a module may assume.
11. **Everything else** [8, 9, 10]: `/search_batch` and quantisation (only
    once per-pixel matching exists -- latency is already 47 ms), class 0 vs
    NULL in the class map, the coverage print mode, the `min_reflectance`
    contrast test, `reflectance_scale` from band metadata rather than
    inference, generic plastic rule floors, the cosmetic `match_sam` EPS.

The OpenCL device-adaptation items in `i.hyper.smac/TODO.md` and
`i.hyper.vegcamopaint/TODO.md` are deliberately outside this ranking: they
are a different subsystem and form a parallel track.

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

The same `library=` serves a second source of reference spectra, now that
`i.hyper.endmembers -i` labels what it extracts: its CSV is one row per
endmember with a `wl_*` column per band and the matched material in
`match_title`, so transposing it into the `wavelength_nm` + `<class>_mean`
shape gives a per-pixel map of the scene's own endmembers rather than of the
field campaign's classes. That is step 2 of the extraction pipeline and is
untouched; the 20 endmembers of the C2RD mosaic are sitting in
`~/grassdata/C2RD_drone_endmembers.csv` waiting for it. No `<class>_sd`
column exists for them, so `library=` must tolerate a missing standard
deviation.

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

- **The served index needs a category filter in practice, not just in the
  API.** Identifying the 20 N-FINDR endmembers of the C2RD drone mosaic
  against the 251,753-spectrum index returned, among plausible vegetation and
  microplastic matches, *Silver*, *Kodak White*, *L-Proline* and *75 vol% Ti
  isopropoxide PDMS* -- laboratory preparations that cannot be in an outdoor
  scene. The cause is library composition, not scoring: `ecosis` and `relab`
  are 96% of the index (210,687 + 31,896 of 251,753), and RELAB's lab mixtures
  match any smooth bright spectrum. `i.hyper.speclookup category=` /
  `i.hyper.endmembers spec_category=` already exist and the server applies
  them, so the work is to find which categories a field scene should be
  restricted to and make that the documented default practice for
  identification, with a measured before/after on the same 20 endmembers.
  A white reference panel is a real possibility in a drone campaign, so
  *Kodak White* is not automatically wrong -- that is exactly what the
  measurement would settle.

- **Nothing cross-checks whether an endmember is pure or merely extreme.**
  `i.hyper.endmembers homogeneity=` now samples a local-variability raster at
  each selected pixel and reports it against the raster's region mean, and
  `i.variance.opencl method=sam_mean|spherical` produces exactly that map for
  a cube. The two have not been run together on the C2RD mosaic, because no
  variance raster exists for `mosaic_final` yet and `i.variance.opencl`
  requires a 3D region matching the map, so the map has to be built at the
  native 0.03 m grid over all 408 bands. Build it, rerun the identification
  with `homogeneity=`, and record which of the 20 endmembers sit in
  heterogeneous neighbourhoods -- those are mixtures or edges whatever their
  library match says.

## 10. A local spectral-lookup service (prototype measured, not yet packaged)

The campaign's lesson was that an *angle against measured spectra* beats a
rule over absorption depths (item 5). The libraries to compare against are
already harvested: `~/grassdata/hyperspeclib` holds **273,002 spectra** in a
Hive-partitioned Parquet tree (ecosis 229k, relab 34k, gsti 7.5k,
usgs_splib07 2.5k), built by the `i.hyper.lib_*` harvesters and queried by
`i.hyper.speclookup`.

A prototype service lives at `~/dev/specserve` (`build_index.py`,
`derive_matrices.py`, `server.py`). Measured on this machine:

| | |
|---|---|
| indexed | 251,753 spectra (21,249 skipped, under 30 % grid coverage) |
| canonical grid | 350-2500 nm at 10 nm, 216 bands |
| resident | 653 MB across 3 memmapped matrices |
| startup | 376 ms |
| query, end to end | **46-69 ms** including JSON and resampling |

### Why there is no tree, and no approximate index

The search is **exact**, a full pass over every candidate, because the cost
is memory bandwidth and not comparisons:

| | time | achieved |
|---|---|---|
| GEMV, 1 query | 6.7 ms | 35.1 GB/s, 17.6 GFLOP/s |
| GEMM, 64 queries | 0.67 ms/query | 177 GFLOP/s |
| GEMM, 256 queries | 0.39 ms/query | 305 GFLOP/s |

A matrix-vector product has about 2 FLOP per byte, so it sits far to the
left of this machine's roofline balance point: the arithmetic is free and
the memory pass is everything. Consequences:

- A space-partitioning tree (KD-tree and friends) reduces *comparisons*,
  which are not the cost, and above roughly 20 dimensions it has to visit
  most of the data anyway. At 216 bands it would pay the same bandwidth,
  lose exactness, and add an index to maintain. Rejected on measurement.
- Batching is the real speedup: Level-2 to Level-3 BLAS reuses each loaded
  cache block across queries and becomes compute-bound, **17x better per
  query at a batch of 256**. A full-map classification must therefore send
  blocks of spectra, not loop single queries. A `/search_batch` endpoint is
  the single highest-value addition.

### Two defects found and fixed in the prototype, worth keeping in mind

- **Do not slice the library per query.** The first version fancy-indexed
  the candidate rows out of the memmap, copying about 200 MB per query:
  740 ms against the 7 ms the arithmetic takes. Replaced by three
  precomputed matrices -- values (zero-filled), their squares, and a 0/1
  support mask -- so the shared-support cosine is four whole-matrix GEMVs
  with no extraction. 10-16x faster.
- **Normalise over the *shared* support, not each spectrum's own.** The
  stored row norm covers a spectrum's full range, which is not its overlap
  with any given query, and reusing it biases the ranking against
  partial-overlap candidates. Fixing it changed the top hit for the dry-soil
  target 072 from a USGS Rangeland soil at 0.060 rad to an EcoSIS spectrum
  at 0.054 rad with 189 of 203 bands shared. This is the same principle as
  item 2: a measurement is only valid over the bands that were measured.

Overlap is enforced rather than assumed -- candidates below `min_overlap`
are rejected, not scored, and every result reports `overlap_bands` against
`query_bands`. On this cube that matters directly: there are no bands at
1329-1482 nm or 1784-1984 nm (item 2), and an angle taken across a gap is
not a measurement.

### Sanity check against surveyed ground

| target | top match | angle |
|---|---|---|
| 072 dry soil excavation | EcoSIS soil, then three USGS *Rangeland* soils | 0.054 rad |
| 097 firethorn shrub | USGS marsh water / vegetation | 0.064 rad |
| 073 transparent polyolefin bag | RELAB organics, USGS **L-Lysine** | 0.089 rad |

All three point the right way, including an organic C--H match for the
polyolefin film -- the target the rule engine could not see at all.

### Still to do

- **Teach the service to carry spectral gaps as no-support -- done.**
  `gaps.py` now holds the segment test, both sides of the comparison use it
  (the query in `Library.search`, the library in `build_index.py`'s
  `GridPlan.inside`), the response carries `query_gaps_nm` and `/health`
  reports `gap_aware`. The measurement that follows is what the fix was
  judged on; the side effect it exposed, and the absolute overlap floor
  added for it, are at the end of this entry. Original diagnosis:
  `Library.search`
  projects a query onto the canonical grid with a single
  `np.interp(grid[inside], q_wl, q_val)` over `[q_wl[0], q_wl[-1]]`, so every
  grid point between the query's first and last wavelength gets a value --
  including the ones the sensor never sampled. A straight line is drawn across
  each water-vapour gap and then compared against candidates that have real
  structure there. This is the same defect as item 2 above, one level up: there
  it was the continuum that had to respect the gaps, here it is the resampling
  onto the shared grid.

  The fix is the same shape too. Split the query into contiguous segments on
  the same relative test `_spectral_segments()` uses (spacing above
  `max(50 nm, 5 x median spacing)`), interpolate within each segment only, and
  leave the gap grid points NaN so they fall out of the support mask -- which
  the scoring already handles, since `qmask`/`q2`/`qm` are exactly the
  machinery for a query that does not cover the whole grid. Nothing in the
  four-GEMV scoring needs to change, and `min_overlap` then measures real
  shared bands instead of counting fabricated ones.

  Measured on one real C2RD drone spectrum (HySpex VS620, 408 bands,
  411.4-2449.7 nm, gaps at 1329.0-1482.6 nm and 1784.5-1984.0 nm) against the
  251,753-spectrum index at `min_overlap=0.55`: **36 of the 203 in-range grid
  points, 17.7%, are fabricated by bridging**. Carrying the gaps instead drops
  the query from 203 bands to 167, and of the top 10 **3 records change** while
  the top 1 holds. Every angle is inflated by the bridging, by 2.0-2.6 degrees
  over the top 10 (best match 20.62 vs 18.61 degrees), because the straight
  line disagrees with whatever the candidate really does in the gap. So the
  bias is conservative -- it makes good matches look worse rather than
  inventing matches -- but it is systematic, it scales with how much of a
  sensor's range is missing, and it reorders near-ties. Any sensor with
  water-vapour gaps is affected: EMIT, EnMAP and PRISMA all are.

  Measured again over all 20 C2RD endmembers once both sides of the index
  were gap-aware, against three servers (old code and index; new code, old
  index; new code, new index): **11 of the 20 top-1 identifications change**,
  mean angle change -0.41 degrees (-1.71 to +1.53), and a mean 4.3 of the
  top 10 records are retained. The library-side half of the fix accounts for
  only 1 of those 11 on this query set -- it excluded 53 spectra and 44 more
  wavelength grids, which had met `min_coverage` only by counting their own
  gaps -- but it is the half that matters for a query with complete coverage,
  where the query-side fix does nothing and only the candidates have gaps.

  **The fraction's denominator was the trap.** `need = min_overlap * qn`, and
  carrying the gaps drops `qn` from 203 to 167, so `min_overlap=0.55` fell
  from a 112-point floor to a 92-point one: the threshold loosened by 18%
  exactly where the query was gappiest. Six of the new top-1 hits rested on
  overlap below the old floor, and restoring it in real bands (112/167 =
  0.67) changed 7 of the 20 again. `Library.search` therefore takes a
  `min_overlap_bands` alongside the fraction and applies the stricter
  (`/health` reports `absolute_overlap_floor`, the response reports
  `overlap_required`); *i.hyper.speclookup* sends both rather than converting
  one into the other, and warns when the server is too old to apply the
  floor. A bare fraction cannot express this requirement, because its
  denominator is a property of the query.

- **Which library is canonical -- settled, and verified end to end
  (2026-11-01).** The preferred path below was taken: openspeclib's
  validated schema is canonical, and the `i.hyper.lib_*` harvest is bridged
  into it by `openspeclib/loaders/hyperspeclib.py`. Verified on the live
  store rather than assumed: all 273,002 records in
  `~/grassdata/openspeclib/catalog.json` carry a `source.filename` of the
  form `dataset_*.parquet`, which only the hyperspeclib tree produces, and
  the per-source counts match `~/grassdata/hyperspeclib/_manifest.json`
  exactly (ecosis 229,368, relab 33,616, gsti 7,561, usgs_splib07 2,457).
  Not one record came from a direct loader. The chain in full:

      $HOME/DBDATA/<download>
        -> i.hyper.lib_{ecosis,gsti,relab,usgs}   (GRASS addon, harvest)
        -> ~/grassdata/hyperspeclib               (Hive Parquet + manifest)
        -> openspeclib loaders/hyperspeclib.py    (validate, licence, schema)
        -> ~/grassdata/openspeclib                (catalog.json + per-source
                                                   parquet, versioned)
        -> specserve build_index.py               (216-point canonical grid)
        -> the served index                       (251,700 of the 273,002)

  The 21,302 dropped at index time are the out-of-grid and unusable-grid
  records, not a second library.

  **Keep it a single entry path.** The reason the verification matters is
  that openspeclib also registers direct loaders for `asu_tes`, `bishop`,
  `ecostress` and `ossl` (`cli.py` `LOADERS`), and those have no
  `i.hyper.lib_*` harvester, so the only way their data can enter is by
  bypassing the GRASS layer. For `ecosis`, `relab` and `usgs` the direct
  loaders are harmless -- `write_source` writes one file per source, so a
  hyperspeclib ingest replaces rather than duplicates them -- but a source
  that enters only through openspeclib is invisible to
  `i.hyper.speclookup library=`, which scans the harvest tree, and has no
  GRASS-side provenance or re-harvest manifest. Rule to hold to: **every
  source gets an `i.hyper.lib_<name>` harvester, and openspeclib's direct
  loader for it is then only a fallback for people without GRASS.**

  Local downloads already in `$HOME/DBDATA` with no harvester, so currently
  unreachable from GRASS (sizes as on disk):

  | Download | Size | openspeclib loader? |
  |---|---|---|
  | `EU_Soil_Spectral_Libraries/OSSL_v1.2` (VNIR only; MIR not downloaded) | 199 MB | yes, `ossl` |
  | `ecospeclib_1770742794330` (ECOSTRESS/ASTER) | 157 MB | yes, `ecostress` |
  | `USDA_ARS_Beltsville` (crop residue, WV3/SWIR) | 56 MB | no |
  | `CannabisSativa` (multispectral, Mendeley) | 23 MB | no |
  | `LUMA_SLUM` (urban materials, IR + SW) | 20 MB | no |
  | `GHISACONUS_2008_001_speclib.csv` + two `GHISACASIA_*.xlsx` (GHISA crops) | 13 MB | no |
  | `CROPS` (DESIS/PRISMA spectral library release) | 6 MB | no |
  | `WaRM_spectral_library` (one CSV per material) | 2 MB | no |
  | `Papers_crop_residues` | 40 KB | no |

  OSSL and ECOSTRESS are the two to do first, because they are the ones a
  direct loader can already pull in behind the harvest's back. OSSL also
  doubles the store's soil coverage, which is currently 5,092 of 273,002
  records (1.9%) against 189,324 vegetation.

- `/search_batch`, for the GEMM speedup above.
- A shared client helper in `i_hyper_lib`, so modules do not each
  reimplement the contract, and so a module degrades gracefully when the
  service is not running instead of failing silently.
- Debian packaging with the index built post-install in the background: a
  systemd **user** service, following the `xerj.service` precedent, since
  the store lives under `$HOME`. Debian policy forbids long work in
  maintainer scripts, so `postinst` registers units and exits while a
  oneshot unit builds the index and reports to MetaMonitor; the server
  answers `/health` with a building status until the matrix exists.
- Quantisation (float16 or int8) is the only lever left on single-query
  latency, since it reduces the bytes read. Validate against the float32
  ranking before adopting.

### Later: how openspeclib gets into GRASS

Not now, but these are the constraints worth recording while they are fresh,
because they decide whether today's structure survives the move.

- **A GRASS addon cannot carry a Python library the way this one is
  carried.** `g.extension` installs a module directory; it does not resolve
  Python dependencies. So openspeclib has to be either an optional
  dependency imported lazily after `gs.parser()` (the pattern every
  `i.hyper.*` module already uses for numpy/pyarrow), with a `gs.fatal`
  naming the pip install when it is missing, or vendored into
  `i_hyper_lib/`. Prefer the optional dependency: it keeps openspeclib's
  own versioning and JSON Schema, which is the whole reason it is canonical.
- **The store does not belong in a mapset.** It is projection-free,
  region-free, read-mostly reference data shared by every project, so
  `$HOME/grassdata/openspeclib` beside the GISDBASE is right and a
  mapset element would be wrong. GRASS has no element type for this, which
  is a real gap: signature files (`struct Signature`) are per-class mean and
  covariance over an imagery group's bands, keyed by semantic label with no
  wavelength axis, and `raster3d` has no semantic-label support at all, so
  neither can hold a library record. Upstream would need a new element, or
  to accept a store outside the mapset.
- **The service is the integration point, not the file format.** Modules
  already talk to `specserve` over HTTP and none of them parses the store,
  so an upstream GRASS could ship the client (item 3 above) and leave the
  store and the server as a separately installed data package. That is also
  what makes the Debian packaging item below the real delivery vehicle.
- **What to propose upstream, in order**: semantic labels (or a wavelength
  axis) for `raster3d`, so a cube can describe itself without `.wl.json`;
  then a spectral-library element or an accepted convention for an
  out-of-mapset store; then the client. The harvesters can stay addons
  indefinitely -- they are downloaders, and upstream has no appetite for
  bundling third-party data.

## Ultimate goal

A fully documented report in
`/home/yann/Documents/LaTex/Report/20261030_ORUS_C2RD_Drone_Flight_Experiment`
with extensive spectral analysis -- literally a spectral forensic analysis.
