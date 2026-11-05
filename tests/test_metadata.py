"""Tests for the processing metadata and the output file name."""

import json

from postprocessing.metadata import (
    build_output_name,
    collect_metadata,
    read_pipeline_versions,
    read_slc_orbits,
)
from postprocessing.run_outputs import find_run_outputs


def test_slcs_are_paired_with_their_orbit_files(run_dir):
    slcs = read_slc_orbits(run_dir / "logs" / "isce")
    assert [slc["date"] for slc in slcs] == ["2022-01-01", "2022-01-13", "2022-01-25"]
    assert all(slc["mission"] == "S1A" and slc["orbit_type"] == "POEORB" for slc in slcs)
    assert slcs[0]["orbit_files"][0].startswith("S1A_OPER_AUX_POEORB")


def test_metadata_and_file_name(run_dir):
    metadata = collect_metadata(find_run_outputs(run_dir), area_name="Wien Nord")
    product = metadata["product"]
    assert product["track"] == 73 and product["swaths"] == ["IW1"]
    assert product["reference_point"]["row"] == 1
    assert metadata["software"]["isce2"] == "2.6.5"
    assert metadata["inputs"]["dem_vertical_datum"] == "egm2008"
    assert metadata["inputs"]["orbit_types"] == {"POEORB": 3}
    assert metadata["settings"]["mintpy_template"]["mintpy.deramp"] == "no"
    json.dumps(metadata, default=str)  # must be serializable
    assert build_output_name(metadata) == "S1_ASC_073_IW1_VV_20220101_20220125_Wien-Nord.parquet"


def test_versions_fall_back_to_logs(run_dir):
    (run_dir / "logs" / "versions.json").unlink()
    assert read_pipeline_versions(run_dir)["dolphin"] == "0.42.8"
