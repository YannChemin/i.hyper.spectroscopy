# Installation — i.hyper.spectroscopy

## Prerequisites

| Requirement | Minimum version | Notes |
|-------------|----------------|-------|
| GRASS GIS | 8.2 | 8.4+ recommended |
| Python | 3.8 | Bundled with GRASS on most platforms |
| NumPy | 1.20 | Usually bundled with GRASS |
| grass.script | (bundled) | GRASS Python scripting library |
| grass.pygrass | (bundled) | Used for fast raster I/O (optional fallback available) |

`libgrass_g3d.so` must be present in `$GISBASE/lib` (standard in all GRASS
installations ≥ 8.0).  The module uses it via ctypes for tile-bulk Z-slice
extraction — no compilation required.

---

## Option A — Install via `g.extension` (recommended)

If the module has been published to the GRASS GIS Addons repository:

```bash
# Start GRASS GIS and open a mapset, then:
g.extension extension=i.hyper.spectroscopy
```

To install from a local source tree (e.g., during development):

```bash
g.extension extension=i.hyper.spectroscopy \
             url=/home/yann/dev/i.hyper.spectroscopy
```

`g.extension` installs the script and HTML manual page and registers the
module so it appears in the GRASS module search.

---

## Option B — Manual installation from source

### 1. Clone or copy the source

```bash
git clone <repository-url> ~/dev/i.hyper.spectroscopy
# or copy the directory to a working location
```

### 2. Set the GRASS build environment

```bash
GRASS_PREFIX=$(grass --config path)
export MODULE_TOPDIR=${GRASS_PREFIX}
```

For developers building against a GRASS source tree:

```bash
export MODULE_TOPDIR=/path/to/grass-source-tree
```

### 3. Compile and install

```bash
cd /home/yann/dev/i.hyper.spectroscopy
make MODULE_TOPDIR=${MODULE_TOPDIR}
make install MODULE_TOPDIR=${MODULE_TOPDIR}
```

This installs:

| File | Destination |
|------|-------------|
| `i.hyper.spectroscopy.py` | `$GISBASE/scripts/i.hyper.spectroscopy` |
| `i.hyper.spectroscopy.html` | `$GISBASE/docs/html/i.hyper.spectroscopy.html` |

### 4. Verify installation

Start GRASS GIS and run:

```bash
i.hyper.spectroscopy --help
```

Expected output begins:

```
Description:
 Per-pixel spectral feature interpretation of hyperspectral imagery (400–2500 nm)

Usage:
 i.hyper.spectroscopy [-nipv] input=name [output=name] [confidence=name]
   [min_band_depth=float] [coordinates=string]
   [--overwrite] [--help] [--verbose] [--quiet] [--ui]
```

---

## Option C — Run directly without installing

For quick testing, execute the script directly from the source directory
inside an active GRASS session:

```bash
# Inside a GRASS terminal session:
python3 /home/yann/dev/i.hyper.spectroscopy/i.hyper.spectroscopy.py \
        input=my_hyperspectral_cube \
        -i
```

---

## Installing the full i.hyper suite

`i.hyper.spectroscopy` reads 3D rasters created by `i.hyper.import` and
works best on atmospherically corrected data.  Recommended installation order
for a complete hyperspectral workflow:

```
i.hyper.import        ← required first (creates 3D raster input format)
i.hyper.atcorr        ← or i.hyper.smac  (atmospheric correction)
i.hyper.continuum     ← continuum removal (sharpens absorption feature depths)
i.hyper.spectroscopy  ← this module
i.hyper.geology       ← geological classification (complementary)
i.hyper.indices       ← supplementary spectral indices
i.hyper.albedo        ← broadband albedo
i.hyper.rgb           ← RGB composites for visualisation
```

---

## Troubleshooting

### `Cannot load libgrass_g3d.so`

The module loads `libgrass_g3d.so` via ctypes for fast Z-slice extraction.
If this fails:

```bash
# Check the library exists:
ls $(grass --config path)/lib/libgrass_g3d.so

# If missing, your GRASS installation may be incomplete.
# On Debian/Ubuntu:
sudo apt-get install grass-dev
```

### `No wavelength metadata found`

`i.hyper.spectroscopy` reads wavelength metadata from two sources (tried in order):

1. **3D raster history** — `r3.info -h` lines formatted as
   `Band N: WAVELENGTH nm, FWHM: F nm`.  Verify with:
   ```bash
   r3.info -h map=my_cube
   ```

2. **Per-band 2D rasters** — individual rasters named `my_cube#N` with
   `r.support` metadata (`wavelength=`, `FWHM=`, `valid=`, `unit=`).
   Verify with:
   ```bash
   r.info -h map="my_cube#1"
   ```

If neither source is present, re-import using `i.hyper.import` or add
metadata manually:

```bash
# Write wavelength metadata for band 1:
r.support map="my_cube#1" \
          description="wavelength=450.0
FWHM=10.0
valid=1
unit=nm"
```

### `Rast3d_extract_z_slice failed`

The fast Z-slice extractor requires the GRASS session's current 3D region to
match the map.  Set it explicitly before running:

```bash
g.region raster_3d=my_cube
i.hyper.spectroscopy input=my_cube output=my_class confidence=my_conf
```

### `grass.pygrass` import error (raster read fallback)

If `grass.pygrass` is unavailable (unusual), the module falls back to
`r.out.ascii` + `numpy.loadtxt`.  This is slower but fully functional.
Install the full GRASS Python bindings to restore fast I/O:

```bash
pip install grass-session   # or reinstall GRASS with Python support
```

### Permissions error during `make install`

If you do not have write access to `$GISBASE`, install into a user-writeable
addon directory:

```bash
export GRASS_ADDON_PATH=$HOME/.grass8/addons
make install MODULE_TOPDIR=${MODULE_TOPDIR} \
             INST_DIR=${GRASS_ADDON_PATH}
```

Then add `GRASS_ADDON_PATH` to your shell environment permanently.

---

## Uninstalling

```bash
# Via g.extension:
g.extension extension=i.hyper.spectroscopy operation=remove

# Manually:
rm -f ${GISBASE}/scripts/i.hyper.spectroscopy
rm -f ${GISBASE}/docs/html/i.hyper.spectroscopy.html
```
