"""Tests for finding the files of a pipeline run."""

import pytest

from postprocessing.run_outputs import RunOutputError, find_final_timeseries, find_run_outputs


def test_final_timeseries_is_the_end_of_the_correction_chain(tmp_path):
    for name in ("timeseries.h5", "timeseries_ERA5.h5", "timeseries_ERA5_demErr.h5",
                 "timeseriesResidual.h5"):
        (tmp_path / name).write_text("")
    assert find_final_timeseries(tmp_path).name == "timeseries_ERA5_demErr.h5"


def test_timeseries_files_outside_one_chain_are_ambiguous(tmp_path):
    for name in ("timeseries_ERA5.h5", "timeseries_GACOS.h5"):
        (tmp_path / name).write_text("")
    with pytest.raises(RunOutputError, match="--timeseries-file"):
        find_final_timeseries(tmp_path)


def test_run_outputs_find_mintpy_and_dolphin_layers(run_dir):
    outputs = find_run_outputs(run_dir)
    assert outputs.shape == (6, 5)
    assert outputs.reference_date.isoformat() == "2022-01-01"
    assert {"lat", "velocity_mintpy", "is_water", "similarity", "is_ps"} <= set(outputs.layers)
    assert "dem_error" not in outputs.layers
    assert len(outputs.crlb_files) == 3


def test_missing_mintpy_directory_is_an_error(tmp_path):
    with pytest.raises(RunOutputError, match="mintpy"):
        find_run_outputs(tmp_path)
