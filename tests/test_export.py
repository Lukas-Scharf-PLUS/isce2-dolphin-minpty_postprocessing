"""Tests for the GeoParquet export."""

import geopandas as gpd
import numpy as np
import pytest

from conftest import VELOCITY_M_PER_YEAR
from postprocessing.export import (
    ExportError,
    export_points,
    los_unit_vector,
    shift_in_meters,
)
from postprocessing.run_outputs import find_run_outputs

WAVELENGTH = 0.05546576


def test_masked_export_writes_only_pixels_in_mintpys_mask(run_dir, tmp_path):
    output = tmp_path / "out.parquet"
    count = export_points(find_run_outputs(run_dir), output, WAVELENGTH)
    points = gpd.read_parquet(output)
    # 6 x 5 grid, first row outside the mask
    assert count == len(points) == 25
    assert points["row"].min() == 1
    assert points.crs.to_string() == "OGC:CRS84"
    assert np.allclose(points.geometry.x, points["lon"])


def test_values_are_converted_to_mm(run_dir, tmp_path):
    output = tmp_path / "out.parquet"
    export_points(find_run_outputs(run_dir), output, WAVELENGTH)
    points = gpd.read_parquet(output)
    assert np.allclose(points["velocity_mintpy"], VELOCITY_M_PER_YEAR * 1000)
    assert np.allclose(points["D20220101"], 0.0)
    expected_last = VELOCITY_M_PER_YEAR * 1000 * 24 / 365.25
    assert np.allclose(points["D20220125"], expected_last, atol=1e-4)
    # CRLB 0.1 rad -> wavelength / (4 pi) * 0.1 in mm
    assert np.allclose(points["D20220113_std"], WAVELENGTH / (4 * np.pi) * 0.1 * 1000)
    # dem_error is not available in this run and therefore left out
    assert "dem_error" not in points.columns


def test_all_pixels_includes_unmasked_but_not_pixels_without_position(run_dir, tmp_path):
    output = tmp_path / "out.parquet"
    export_points(find_run_outputs(run_dir), output, WAVELENGTH, all_pixels=True)
    points = gpd.read_parquet(output)
    assert len(points) == 29  # 30 pixels, pixel (0, 0) has no position
    water = points[(points["row"] == 5) & (points["col"] == 4)]
    assert water["is_water"].item() == 1
    assert water["water_fraction"].item() == pytest.approx(0.9)


def test_excluded_groups_are_left_out(run_dir, tmp_path):
    output = tmp_path / "out.parquet"
    export_points(find_run_outputs(run_dir), output, WAVELENGTH,
                  excluded=["displacement_std", "model", "phase_linking"])
    columns = set(gpd.read_parquet(output).columns)
    assert "D20220113_std" not in columns and "mean_velocity" not in columns
    assert "similarity" not in columns and "D20220113" in columns


def test_existing_output_needs_overwrite(run_dir, tmp_path):
    output = tmp_path / "out.parquet"
    output.write_text("")
    with pytest.raises(ExportError, match="overwrite"):
        export_points(find_run_outputs(run_dir), output, WAVELENGTH)


def test_los_unit_vector_points_west_and_up_for_ascending():
    east, north, up = los_unit_vector(np.array([35.0]), np.array([102.0]))
    assert east[0] < 0 and up[0] > 0
    assert east[0] ** 2 + north[0] ** 2 + up[0] ** 2 == pytest.approx(1.0)


def _add_subpixel_output(run_dir, shape=(6, 5)):
    """Write the pipeline's subpixel/ GeoTIFFs: one PS at (2, 3), NaN elsewhere.

    The PS sits 1e-5 deg north and 2e-5 deg east of its pixel center.
    """
    from conftest import write_tif
    import h5py

    with h5py.File(run_dir / "mintpy" / "inputs" / "geometryRadar.h5", "r") as geometry:
        lat = geometry["latitude"][()].astype(np.float64)
        lon = geometry["longitude"][()].astype(np.float64)
    layers = {name: np.full(shape, np.nan) for name in (
        "latitude", "longitude", "height", "range_offset", "azimuth_offset",
        "range_offset_std", "azimuth_offset_std", "valid_dates", "on_peak")}
    values = {"latitude": lat[2, 3] + 1e-5, "longitude": lon[2, 3] + 2e-5, "height": 201.0,
              "range_offset": 0.3, "azimuth_offset": -0.4, "range_offset_std": 0.05,
              "azimuth_offset_std": 0.06, "valid_dates": 3, "on_peak": 1}
    for name, value in values.items():
        layers[name][2, 3] = value
        dtype = np.float64 if name in ("latitude", "longitude") else np.float32
        write_tif(run_dir / "subpixel" / f"{name}.tif", layers[name].astype(dtype))


def test_subpixel_columns_are_written_when_the_step_ran(run_dir, tmp_path):
    _add_subpixel_output(run_dir)
    output = tmp_path / "out.parquet"
    export_points(find_run_outputs(run_dir), output, WAVELENGTH)
    points = gpd.read_parquet(output)

    ps = points[(points["row"] == 2) & (points["col"] == 3)].iloc[0]
    assert ps["lat_subpixel"] == pytest.approx(ps["lat"] + 1e-5, abs=1e-6)
    assert ps["shift_north"] == pytest.approx(1e-5 * 111_200, rel=0.01)
    assert ps["shift_east"] == pytest.approx(2e-5 * 74_300, rel=0.01)
    assert ps["subpixel_on_peak"] == 1 and ps["subpixel_valid_dates"] == 3
    assert ps["azimuth_offset"] == pytest.approx(-0.4)
    # Original position columns stay the pixel center; the point geometry too.
    assert ps.geometry.y == pytest.approx(ps["lat"])
    # Pixels without PS have empty sub-pixel columns.
    others = points[~((points["row"] == 2) & (points["col"] == 3))]
    assert others["lat_subpixel"].isna().all() and others["shift_east"].isna().all()


def test_no_subpixel_output_leaves_the_columns_out_without_warning(run_dir, tmp_path, caplog):
    output = tmp_path / "out.parquet"
    export_points(find_run_outputs(run_dir), output, WAVELENGTH)
    columns = set(gpd.read_parquet(output).columns)
    assert not columns & {"lat_subpixel", "shift_east", "range_offset", "subpixel_on_peak"}
    assert not any("lat_subpixel" in record.message for record in caplog.records
                   if record.levelname == "WARNING")


def test_subpixel_group_can_be_excluded(run_dir, tmp_path):
    _add_subpixel_output(run_dir)
    output = tmp_path / "out.parquet"
    export_points(find_run_outputs(run_dir), output, WAVELENGTH, excluded=["subpixel"])
    assert "shift_east" not in gpd.read_parquet(output).columns


def test_dem_error_is_written_when_mintpy_estimated_it(run_dir, tmp_path):
    from conftest import write_h5
    write_h5(run_dir / "mintpy" / "demErr.h5", {"dem": np.full((6, 5), 2.5, dtype=np.float32)})
    output = tmp_path / "out.parquet"
    export_points(find_run_outputs(run_dir), output, WAVELENGTH)
    assert np.allclose(gpd.read_parquet(output)["dem_error"], 2.5)


def test_shift_in_meters_matches_known_distances():
    east, north = shift_in_meters(np.array([48.2]), np.array([16.3]),
                                  np.array([48.2 + 1e-5]), np.array([16.3 + 1e-5]))
    # At 48.2 deg: 1 deg latitude ~ 111.2 km, 1 deg longitude ~ 74.3 km.
    assert north[0] == pytest.approx(1.1120, abs=0.001)
    assert east[0] == pytest.approx(0.7429, abs=0.001)
