#!/usr/bin/env python
"""Command-line entry point: export an isce2-dolphin-mintpy run to one GeoParquet point file.

The file is named after its content, e.g.
S1_ASC_073_IW1_VV_20220101_20220125_vienna.parquet, and contains the processing
metadata (also written next to it as <name>.metadata.json).

Usage:
  python run_postprocessing.py --run-dir <pipeline work_dir> --area-name vienna
  python run_postprocessing.py --run-dir <dir> --output-dir <folder>
  python run_postprocessing.py --run-dir <dir> --all-pixels
  python run_postprocessing.py --run-dir <dir> --exclude displacement_std
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import h5py

from postprocessing.attributes import FIXED_ATTRIBUTES, GROUPS
from postprocessing.export import DEFAULT_BLOCK_ROWS, ExportError, export_points
from postprocessing.metadata import (
    METADATA_KEY,
    build_output_name,
    clean_area_name,
    collect_metadata,
    write_sidecar,
)
from postprocessing.run_outputs import RunOutputError, find_run_outputs

logger = logging.getLogger("postprocessing")

# Sentinel-1 C-band, used if the time series has no WAVELENGTH attribute.
SENTINEL1_WAVELENGTH_M = 0.055465764662349676


def create_parser() -> argparse.ArgumentParser:
    """Create the command-line parser."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--run-dir", type=Path, required=True,
        help="The pipeline's work directory (contains mintpy/ and dolphin/).")
    parser.add_argument(
        "--timeseries-file", type=Path,
        help="MintPy time series to export. Default: the end of MintPy's correction chain.")
    output = parser.add_mutually_exclusive_group()
    output.add_argument(
        "--output-dir", type=Path,
        help="Folder for the output file, which is named after its content. Default: --run-dir.")
    output.add_argument("--output", type=Path, help="Output file, to choose the name yourself.")
    parser.add_argument(
        "--area-name", help="Name of the area, e.g. vienna. Added to the file name and metadata.")
    parser.add_argument(
        "--exclude", nargs="+", default=[], metavar="NAME",
        help=f"Groups or columns to leave out. Groups: {', '.join(GROUPS)}. "
        f"Columns: {', '.join(a.name for a in FIXED_ATTRIBUTES)}.")
    parser.add_argument(
        "--all-pixels", action="store_true",
        help="Write every pixel with a position, not only those in MintPy's mask "
        "(maskTempCoh.h5). The mask is still written as a column.")
    parser.add_argument(
        "--wavelength", type=float,
        help="Radar wavelength in meters. Default: from the time series (Sentinel-1).")
    parser.add_argument(
        "--block-rows", type=int, default=DEFAULT_BLOCK_ROWS,
        help="Radar rows processed at a time. Smaller uses less memory.")
    parser.add_argument("--overwrite", action="store_true", help="Replace an existing output file.")
    return parser


def _wavelength(args: argparse.Namespace, timeseries_file: Path) -> float:
    """Return --wavelength, else the time series' WAVELENGTH, else Sentinel-1."""
    if args.wavelength:
        return args.wavelength
    with h5py.File(timeseries_file, "r") as h5_file:
        value = h5_file.attrs.get("WAVELENGTH")
    return float(value) if value is not None else SENTINEL1_WAVELENGTH_M


def main() -> int:
    """Parse the command line and write the point file.

    Returns
    -------
    int
        Process exit code: 0 on success, 1 if the export failed.
    """
    args = create_parser().parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s")
    run_dir = args.run_dir.expanduser().resolve()
    timeseries_file = args.timeseries_file.expanduser().resolve() if args.timeseries_file else None
    # An empty name (e.g. an empty Argo parameter) means no area name.
    area_name = args.area_name or None

    # Expected problems get a clear message instead of a long stack trace.
    try:
        if area_name is not None:
            clean_area_name(area_name)  # fail before the slow export
        outputs = find_run_outputs(run_dir, timeseries_file)
        wavelength_m = _wavelength(args, outputs.timeseries_file)
        metadata = collect_metadata(outputs, area_name=area_name, wavelength_m=wavelength_m,
                                    excluded=args.exclude, all_pixels=args.all_pixels)
        output_file = (args.output.expanduser().resolve() if args.output else
                       (args.output_dir or run_dir).expanduser().resolve()
                       / build_output_name(metadata))
        export_points(
            outputs, output_file, wavelength_m=wavelength_m, excluded=args.exclude,
            all_pixels=args.all_pixels, block_rows=args.block_rows, overwrite=args.overwrite,
            file_metadata={METADATA_KEY: json.dumps(metadata, default=str)},
        )
        logger.info("Wrote the metadata also to %s", write_sidecar(metadata, output_file))
    except (RunOutputError, ExportError, ValueError) as error:
        logger.error("%s", error)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
