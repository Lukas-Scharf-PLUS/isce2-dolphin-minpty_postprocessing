"""Tests for the GeoParquet export."""

import geopandas as gpd
import numpy as np
import pytest

from conftest import VELOCITY_M_PER_YEAR
from postprocessing.export import ExportError, export_points, los_unit_vector
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
