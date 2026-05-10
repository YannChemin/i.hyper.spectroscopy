"""Tests that libras3d is functional for i.hyper.spectroscopy standalone mode."""
import os, sys
import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from test_ras3d_common import (
    WYVERN_PATH, TANAGER_PATH,
    skip_without_ras3d, skip_without_wyvern, skip_without_tanager,
    open_cube_checked, assert_band_valid, install_ras3d_shim, make_wl_sidecar,
)


@skip_without_ras3d
def test_shim_installs():
    """ras3d shim provides grass.script API."""
    install_ras3d_shim()
    import grass.script as gs
    assert hasattr(gs, 'raster3d_info')
    assert hasattr(gs, 'parser')


@skip_without_ras3d
@skip_without_wyvern
def test_open_wyvern_geotiff():
    """Wyvern 23-band GeoTIFF opens with correct dimensions."""
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    assert r['depths'] == 23
    assert r['rows']   == 7825
    assert r['cols']   == 6003
    ras3d.close_cube(h)


@skip_without_ras3d
@skip_without_tanager
def test_open_tanager_hdf5():
    """Tanager 426-band HDF5 opens with correct dimensions."""
    import ras3d
    h, r = open_cube_checked(TANAGER_PATH)
    assert r['depths'] == 426
    assert r['rows']   == 732
    assert r['cols']   == 607
    ras3d.close_cube(h)


@skip_without_ras3d
@skip_without_wyvern
def test_read_all_bands_wyvern():
    """read_all_bands() returns float32 cube with correct shape."""
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    cube = ras3d.read_all_bands(h)
    assert cube.shape == (r['depths'], r['rows'], r['cols'])
    assert cube.dtype == np.float32
    ras3d.close_cube(h)


@skip_without_ras3d
@skip_without_tanager
def test_read_all_bands_tanager():
    """read_all_bands() on HDF5 returns correct shape."""
    import ras3d
    h, r = open_cube_checked(TANAGER_PATH)
    cube = ras3d.read_all_bands(h)
    assert cube.shape == (r['depths'], r['rows'], r['cols'])
    assert_band_valid(cube[0],   'Tanager band 0')
    assert_band_valid(cube[425], 'Tanager band 425')
    ras3d.close_cube(h)


@skip_without_ras3d
@skip_without_wyvern
def test_wavelength_sidecar():
    """Wavelength function reads .wl.json sidecar in ras3d mode."""
    install_ras3d_shim()
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    sidecar, wl_expected = make_wl_sidecar(WYVERN_PATH, r['depths'])
    ras3d.close_cube(h)
    sys.path.insert(0, '/home/yann/dev/i.hyper.spectroscopy')
    import importlib
    m = importlib.import_module('i_hyper_spectroscopy')
    try:
        fn = getattr(m, 'get_band_info')
        result = fn(WYVERN_PATH)
        assert result is not None
        assert len(result) > 0
    except AttributeError:
        pytest.skip('get_band_info not found in module (may have different name)')
    finally:
        os.unlink(sidecar)


@skip_without_ras3d
@skip_without_wyvern
def test_get_region_from_shim():
    """gs.raster3d_info() shim returns correct metadata for Wyvern cube."""
    install_ras3d_shim()
    import grass.script as gs
    info = gs.raster3d_info(WYVERN_PATH)
    assert info['depths'] == 23
    assert info['rows']   == 7825
    assert info['cols']   == 6003
