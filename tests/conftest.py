"""Shared test helpers: a tiny isce2-dolphin-mintpy run directory with synthetic data."""

import datetime
import json

import h5py
import numpy as np
import pytest
import rasterio
import yaml

SHAPE = (6, 5)
DATES = [datetime.date(2022, 1, 1), datetime.date(2022, 1, 13), datetime.date(2022, 1, 25)]
VELOCITY_M_PER_YEAR = -0.01
SAFES = [
    f"S1A_IW_SLC__1SDV_{date:%Y%m%d}T164305_{date:%Y%m%d}T164333_041445_04ED7B_5516.SAFE"
    for date in DATES
]
ORBITS = [
    f"S1A_OPER_AUX_POEORB_OPOD_20220202T081558_V{date:%Y%m%d}T000000_{date:%Y%m%d}T235959.EOF"
    for date in DATES
]
ATTRIBUTES = {
    "FILE_TYPE": "timeseries", "UNIT": "m", "REF_DATE": "20220101", "REF_Y": "1", "REF_X": "2",
    "WAVELENGTH": "0.05546576", "ORBIT_DIRECTION": "ASCENDING", "trackNumber": "73",
    "POLARIZATION": "VV", "swathNumber": "1", "RANGE_PIXEL_SIZE": "9.3",
    "AZIMUTH_PIXEL_SIZE": "15.6",
}


def write_h5(path, datasets, attributes=None):
    """Write an HDF5 file with datasets and file attributes, creating its folder."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h5_file:
        for name, values in datasets.items():
            h5_file[name] = values
        for key, value in (attributes or {}).items():
            h5_file.attrs[key] = value
    return path


def write_tif(path, values):
    """Write a one-band GeoTIFF without georeferencing (radar grid)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", driver="GTiff", height=values.shape[0], width=values.shape[1],
                       count=1, dtype=values.dtype) as raster:
        raster.write(values, 1)
    return path


@pytest.fixture
def run_dir(tmp_path):
    """A pipeline run: MintPy files, dolphin rasters, ISCE2 logs and configs.

    Displacement is a linear motion of VELOCITY_M_PER_YEAR. In maskTempCoh, the
    first row is False. Pixel (0, 0) has no position (lat = lon = 0).
    """
    rows, cols = np.mgrid[0:SHAPE[0], 0:SHAPE[1]]
    years = np.array([(date - DATES[0]).days / 365.25 for date in DATES])
    timeseries = (years[:, None, None] * VELOCITY_M_PER_YEAR * np.ones(SHAPE)).astype(np.float32)
    latitude = (48.2 - rows * 0.0001).astype(np.float32)
    longitude = (16.3 + cols * 0.0002).astype(np.float32)
    latitude[0, 0] = longitude[0, 0] = 0.0
    mask = np.ones(SHAPE, dtype=bool)
    mask[0] = False
    water_fraction = np.zeros(SHAPE, dtype=np.float32)
    water_fraction[-1, -1] = 0.9

    mintpy = tmp_path / "mintpy"
    write_h5(mintpy / "timeseries.h5", {
        "timeseries": timeseries,
        "date": np.array([f"{date:%Y%m%d}" for date in DATES], dtype="S8"),
        "bperp": np.zeros(len(DATES), dtype=np.float32),
    }, ATTRIBUTES)
    write_h5(mintpy / "inputs" / "geometryRadar.h5", {
        "latitude": latitude, "longitude": longitude,
        "height": np.full(SHAPE, 200.0, dtype=np.float32),
        "incidenceAngle": np.full(SHAPE, 35.0, dtype=np.float32),
        "azimuthAngle": np.full(SHAPE, 102.0, dtype=np.float32),
        "shadowMask": np.zeros(SHAPE, dtype=bool),
        "waterFraction": water_fraction, "waterMask": water_fraction <= 0.8,
    })
    write_h5(mintpy / "velocity.h5", {
        "velocity": np.full(SHAPE, VELOCITY_M_PER_YEAR, dtype=np.float32),
        "velocityStd": np.full(SHAPE, 0.001, dtype=np.float32),
    })
    write_h5(mintpy / "temporalCoherence.h5",
             {"temporalCoherence": np.full(SHAPE, 0.9, dtype=np.float32)})
    write_h5(mintpy / "maskTempCoh.h5", {"mask": mask})
    write_h5(mintpy / "dolphinTemporalCoherence.h5",
             {"temporalCoherence": np.full(SHAPE, 0.8, dtype=np.float32)})
    (mintpy / "isce_dolphin_mintpy.cfg").write_text(
        "mintpy.troposphericDelay.method = no\nmintpy.deramp = no\n")

    interferograms = tmp_path / "dolphin" / "interferograms"
    write_tif(interferograms / "similarity_20220101_20220125.tif", np.full(SHAPE, 0.7, np.float32))
    write_tif(interferograms / "ps_mask_looked.tif", np.zeros(SHAPE, np.uint8))
    for date in DATES:
        write_tif(interferograms / f"crlb_{date:%Y%m%d}.tif", np.full(SHAPE, 0.1, np.float32))
    (tmp_path / "dolphin" / "dolphin_config.yaml").write_text(
        yaml.safe_dump({"output_options": {"strides": {"x": 4, "y": 1}}}))
    (tmp_path / "dolphin" / "dolphin.log").write_text(
        json.dumps({"message": "Current running dolphin version: 0.42.8"}) + "\n")

    isce_logs = tmp_path / "logs" / "isce"
    isce_logs.mkdir(parents=True)
    for index, (safe, orbit) in enumerate(zip(SAFES, ORBITS)):
        (isce_logs / f"run_02_unpack_secondary_slc_{index:03d}.log").write_text(
            f"['--dirname', '/data/slc/{safe}', '--swaths', '1']\n"
            f"Extracting orbit from Orbit File:  /data/orbits/{orbit}\n")
    (tmp_path / "logs" / "versions.json").write_text(json.dumps(
        {"isce2": "2.6.5", "dolphin": "0.42.8", "mintpy": "1.6.4", "isce_dolphin_mintpy": "abc"}))
    (tmp_path / "pipeline_config.yaml").write_text(yaml.safe_dump({
        "paths": {"dem_file": "/data/dem.tif", "dem_vertical_datum": "egm2008"},
        "water_mask": {"file": "/data/watermask.tif"},
    }))
    (tmp_path / "isce" / "reference").mkdir(parents=True)
    (tmp_path / "isce" / "reference" / "IW1.xml").write_text("<x/>")
    return tmp_path
