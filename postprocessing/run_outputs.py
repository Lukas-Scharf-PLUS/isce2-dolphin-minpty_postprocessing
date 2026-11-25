"""Find the files of an isce2-dolphin-mintpy run that become point attributes.

A run directory (the pipeline's ``paths.work_dir``) looks like this::

    <run_dir>/
    ├── pipeline_config.yaml          copy of the pipeline config
    ├── logs/versions.json            software versions (newer pipeline versions)
    ├── isce/                         ISCE2 topsStack (SAFE_files.txt, reference/, baselines/)
    ├── dolphin/                      dolphin work directory (interferograms/, dolphin.log)
    ├── subpixel/                     optional: sub-pixel PS positions (pipeline step "subpixel")
    └── mintpy/                       MintPy work directory, radar coordinates
        ├── inputs/geometryRadar.h5   lat, lon, height, angles, shadow and water layers
        ├── timeseries*.h5            time series; the corrected ones add suffixes
        ├── velocity.h5, temporalCoherence.h5, maskTempCoh.h5, ...

Everything is in **radar coordinates** on dolphin's output grid: MintPy's files
and dolphin's rasters have the same shape, pixel for pixel. Every pixel gets its
exact latitude and longitude from the ISCE2 geometry, so no resampling is needed.

**Final time series.** MintPy writes one file per correction step, each name
extending the previous one, e.g. ``timeseries.h5`` -> ``timeseries_ERA5.h5`` ->
``timeseries_ERA5_demErr.h5``. The final one is the longest name of that chain.

Only the geometry and the time series are required; every other file is optional.
"""

import datetime
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import h5py

logger = logging.getLogger(__name__)

DATE_FORMAT = "%Y%m%d"
DATE_PAIR_PATTERN = re.compile(r"\d{8}_\d{8}")
CRLB_PATTERN = re.compile(r"crlb_(?P<date>\d{8})")

MINTPY_DIRECTORY = "mintpy"
DOLPHIN_DIRECTORY = "dolphin"
INTERFEROGRAM_DIRECTORY = "interferograms"
SUBPIXEL_DIRECTORY = "subpixel"
# Column name -> GeoTIFF written by the pipeline's subpixel step (on dolphin's output grid).
SUBPIXEL_FILES = {
    "lon_subpixel": "longitude.tif",
    "lat_subpixel": "latitude.tif",
    "height_subpixel": "height.tif",
    "range_offset": "range_offset.tif",
    "azimuth_offset": "azimuth_offset.tif",
    "range_offset_std": "range_offset_std.tif",
    "azimuth_offset_std": "azimuth_offset_std.tif",
    "subpixel_valid_dates": "valid_dates.tif",
    "subpixel_on_peak": "on_peak.tif",
}
TIMESERIES_PREFIX = "timeseries"
# MintPy files starting with "timeseries" that are not displacement time series.
NOT_TIMESERIES_PREFIXES = ("timeseriesResidual", "timeseriesStd")

# dolphin ministack products: file prefix, and the tag of the combined file.
SHP_COUNTS = ("shp_counts", "average")
SIMILARITY = ("similarity", "full")


class RunOutputError(RuntimeError):
    """The run output is missing or cannot be interpreted."""


@dataclass(frozen=True)
class Layer:
    """One 2D raster layer: a dataset in an HDF5 file, or band 1 of a GeoTIFF.

    Attributes
    ----------
    path : Path
        The file.
    dataset : str or None
        HDF5 dataset name; None for a GeoTIFF.
    """

    path: Path
    dataset: str | None = None


@dataclass(frozen=True)
class RunOutputs:
    """All files of one run that the export can use.

    Attributes
    ----------
    run_dir : Path
        The pipeline work directory.
    timeseries_file : Path
        MintPy's final (fully corrected) time series.
    dates : list[datetime.date]
        Dates of the time series.
    reference_date : datetime.date
        MintPy's reference date (zero displacement).
    shape : tuple[int, int]
        (rows, columns) of the radar grid.
    layers : dict[str, Layer]
        Optional single-layer sources, by column name (see export.py).
    crlb_files : dict[datetime.date, Path]
        dolphin's per-date CRLB rasters; may be empty.
    full_range_pixel_size_m : float or None
        Slant-range spacing of the full-resolution SLC samples in meters, i.e.
        the unit of the sub-pixel ``range_offset``. None if unknown.
    """

    run_dir: Path
    timeseries_file: Path
    dates: list[datetime.date]
    reference_date: datetime.date
    shape: tuple[int, int]
    layers: dict[str, Layer] = field(default_factory=dict)
    crlb_files: dict[datetime.date, Path] = field(default_factory=dict)
    full_range_pixel_size_m: float | None = None

    @property
    def mintpy_dir(self) -> Path:
        """MintPy's work directory."""
        return self.run_dir / MINTPY_DIRECTORY

    @property
    def dolphin_dir(self) -> Path:
        """dolphin's work directory."""
        return self.run_dir / DOLPHIN_DIRECTORY


def parse_date(text: str) -> datetime.date:
    """Parse a ``YYYYMMDD`` date."""
    return datetime.datetime.strptime(text, DATE_FORMAT).date()


def find_final_timeseries(mintpy_dir: Path) -> Path:
    """Return MintPy's final time series: the end of the correction chain.

    Parameters
    ----------
    mintpy_dir : Path
        MintPy's work directory.

    Returns
    -------
    Path
        The time series file whose name all other time series names lead to.

    Raises
    ------
    RunOutputError
        If there is no time series, or the files do not form one chain (then
        it is unclear which is final; pass ``--timeseries-file``).
    """
    candidates = [
        path for path in mintpy_dir.glob(f"{TIMESERIES_PREFIX}*.h5")
        if not path.name.startswith(NOT_TIMESERIES_PREFIXES)
    ]
    if not candidates:
        raise RunOutputError(f"No {TIMESERIES_PREFIX}*.h5 in {mintpy_dir}.")
    final = max(candidates, key=lambda path: len(path.stem))
    not_in_chain = [path.name for path in candidates if not final.stem.startswith(path.stem)]
    if not_in_chain:
        raise RunOutputError(
            f"The time series files in {mintpy_dir} are not one correction chain "
            f"(final candidate {final.name}, others {not_in_chain}). "
            "Choose one with --timeseries-file."
        )
    return final


def read_full_range_pixel_size(timeseries_file: Path) -> float | None:
    """Return the full-resolution slant-range sample spacing in meters, or None.

    MintPy's ``RANGE_PIXEL_SIZE`` is the spacing of the output grid. The pipeline
    sets it to the full-resolution spacing times dolphin's range stride, which it
    stores as ``RLOOKS`` (1 if missing), so the full-resolution spacing is their
    ratio.
    """
    with h5py.File(timeseries_file, "r") as h5_file:
        pixel_size = h5_file.attrs.get("RANGE_PIXEL_SIZE")
        range_looks = h5_file.attrs.get("RLOOKS", 1)
    if pixel_size is None:
        return None
    return float(_text(pixel_size)) / float(_text(range_looks))


def read_timeseries_dates(timeseries_file: Path) -> tuple[list[datetime.date], datetime.date, tuple]:
    """Read the dates, reference date and grid shape of a MintPy time series.

    Raises
    ------
    RunOutputError
        If the file is not a MintPy time series.
    """
    with h5py.File(timeseries_file, "r") as h5_file:
        if "timeseries" not in h5_file or "date" not in h5_file:
            raise RunOutputError(f"{timeseries_file} is not a MintPy time series.")
        dates = [parse_date(value.decode()) for value in h5_file["date"][()]]
        shape = h5_file["timeseries"].shape[1:]
        reference = h5_file.attrs.get("REF_DATE")
    reference_date = parse_date(_text(reference)) if reference is not None else dates[0]
    return dates, reference_date, tuple(int(size) for size in shape)


def _text(value) -> str:
    """Return an HDF5 attribute (bytes or str) as str."""
    return value.decode() if isinstance(value, bytes) else str(value)


def find_ministack_product(interferogram_dir: Path, prefix: str, combined_tag: str) -> Path | None:
    """Find a dolphin phase linking product, combined over all ministacks.

    Returns ``<prefix>_<combined_tag>_<span>.tif`` if it exists, otherwise the
    only ``<prefix>_<span>.tif``, otherwise None. Several per-ministack files
    without a combined one give a warning and None, because it is unclear which
    to use (same rule as dolphin_postprocessing).
    """
    combined = sorted(interferogram_dir.glob(f"{prefix}_{combined_tag}_*.tif"))
    if combined:
        return combined[-1]
    single = [
        path for path in sorted(interferogram_dir.glob(f"{prefix}_*.tif"))
        if DATE_PAIR_PATTERN.fullmatch(path.stem.removeprefix(f"{prefix}_"))
    ]
    if len(single) > 1:
        logger.warning("Several %s rasters but no combined one in %s; leaving them out.",
                       prefix, interferogram_dir)
        return None
    return single[0] if single else None


def find_crlb_files(interferogram_dir: Path) -> dict[datetime.date, Path]:
    """Find dolphin's per-date CRLB rasters ``crlb_<date>.tif``, sorted by date."""
    crlb_files = {}
    for path in interferogram_dir.glob("crlb_*.tif"):
        match = CRLB_PATTERN.fullmatch(path.stem)
        if match:
            crlb_files[parse_date(match["date"])] = path
    return dict(sorted(crlb_files.items()))


def _h5_layer(path: Path, dataset: str) -> Layer | None:
    """Return a Layer if the HDF5 file exists and has the dataset, else None."""
    if not path.is_file():
        return None
    with h5py.File(path, "r") as h5_file:
        if dataset not in h5_file:
            return None
    return Layer(path, dataset)


def _first_dataset_layer(path: Path) -> Layer | None:
    """Return a Layer for the only dataset of a single-layer MintPy file, or None."""
    if not path.is_file():
        return None
    with h5py.File(path, "r") as h5_file:
        names = [name for name in h5_file if isinstance(h5_file[name], h5py.Dataset)]
    return Layer(path, names[0]) if len(names) == 1 else None


def _tif_layer(path: Path | None) -> Layer | None:
    """Return a Layer for an existing GeoTIFF, else None."""
    return Layer(path) if path is not None and path.is_file() else None


def find_run_outputs(run_dir: Path, timeseries_file: Path | None = None) -> RunOutputs:
    """Find all usable files of a pipeline run.

    Parameters
    ----------
    run_dir : Path
        The pipeline work directory (contains ``mintpy/`` and ``dolphin/``).
    timeseries_file : Path, optional
        Use this time series instead of the end of MintPy's correction chain.

    Returns
    -------
    RunOutputs

    Raises
    ------
    RunOutputError
        If the MintPy directory, its geometry or the time series is missing.
    """
    mintpy_dir = run_dir / MINTPY_DIRECTORY
    if not mintpy_dir.is_dir():
        raise RunOutputError(f"No {MINTPY_DIRECTORY}/ directory in {run_dir}.")
    geometry = mintpy_dir / "inputs" / "geometryRadar.h5"
    if not geometry.is_file():
        raise RunOutputError(f"MintPy geometry missing: {geometry}")

    timeseries_file = timeseries_file or find_final_timeseries(mintpy_dir)
    dates, reference_date, shape = read_timeseries_dates(timeseries_file)
    logger.info("Time series: %s (%d dates, reference date %s)",
                timeseries_file.name, len(dates), reference_date)

    interferogram_dir = run_dir / DOLPHIN_DIRECTORY / INTERFEROGRAM_DIRECTORY
    candidates = {
        "lat": _h5_layer(geometry, "latitude"),
        "lon": _h5_layer(geometry, "longitude"),
        "height": _h5_layer(geometry, "height"),
        "incidence_angle": _h5_layer(geometry, "incidenceAngle"),
        "azimuth_angle": _h5_layer(geometry, "azimuthAngle"),
        "shadow_layover": _h5_layer(geometry, "shadowMask"),
        "water_fraction": _h5_layer(geometry, "waterFraction"),
        # waterMask is True for land; export.py turns it into is_water.
        "is_water": _h5_layer(geometry, "waterMask"),
        "temporal_coherence": _h5_layer(mintpy_dir / "temporalCoherence.h5", "temporalCoherence"),
        "avg_spatial_coherence": _h5_layer(mintpy_dir / "avgSpatialCoh.h5", "coherence"),
        "mask_temporal_coherence": _h5_layer(mintpy_dir / "maskTempCoh.h5", "mask"),
        "conncomp_ok": _h5_layer(mintpy_dir / "maskConnComp.h5", "mask"),
        "num_inverted_ifgrams": _h5_layer(mintpy_dir / "numInvIfgram.h5", "mask"),
        "phase_linking_coherence": _h5_layer(mintpy_dir / "dolphinTemporalCoherence.h5",
                                             "temporalCoherence"),
        "velocity_mintpy": _h5_layer(mintpy_dir / "velocity.h5", "velocity"),
        "velocity_mintpy_std": _h5_layer(mintpy_dir / "velocity.h5", "velocityStd"),
        "dem_error": _first_dataset_layer(mintpy_dir / "demErr.h5"),
        "similarity": _tif_layer(find_ministack_product(interferogram_dir, *SIMILARITY)),
        "shp_count": _tif_layer(find_ministack_product(interferogram_dir, *SHP_COUNTS)),
        "is_ps": _tif_layer(interferogram_dir / "ps_mask_looked.tif"),
        "amp_dispersion": _tif_layer(interferogram_dir / "amp_dispersion_looked.tif"),
    }
    subpixel_dir = run_dir / SUBPIXEL_DIRECTORY
    candidates |= {name: _tif_layer(subpixel_dir / file_name)
                   for name, file_name in SUBPIXEL_FILES.items()}
    for required in ("lat", "lon"):
        if candidates[required] is None:
            raise RunOutputError(f"{geometry} has no latitude/longitude.")
    layers = {name: layer for name, layer in candidates.items() if layer is not None}

    return RunOutputs(
        run_dir=run_dir,
        timeseries_file=timeseries_file,
        dates=dates,
        reference_date=reference_date,
        shape=shape,
        layers=layers,
        crlb_files=find_crlb_files(interferogram_dir) if interferogram_dir.is_dir() else {},
        full_range_pixel_size_m=read_full_range_pixel_size(timeseries_file),
    )
