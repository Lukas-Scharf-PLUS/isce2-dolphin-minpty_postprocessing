"""Collect processing metadata for the point file, and build its file name.

The metadata says what the file contains and how it was made: sensor, track,
dates, the SLC and orbit files used, all software versions and the settings of
every stage. It is stored inside the Parquet file (and as a readable JSON file
next to it), so it cannot get separated from the data.

Sources in the run directory, all optional except the time series:

- MintPy's time series attributes: track, direction, polarization, reference
  date and point,
- ``pipeline_config.yaml``: the complete pipeline config (DEM, datum, AOI, ...),
- ``logs/versions.json``: ISCE2, dolphin, MintPy and pipeline versions
  (fallback: the dolphin and MintPy logs),
- ``logs/isce/run_0[12]_*.log``: each SAFE with the orbit file ISCE2 used,
- ``dolphin/dolphin_config.yaml`` and MintPy's template: their settings.

If a source is missing, its fields are None and a warning is logged.
"""

import datetime
import importlib.metadata
import json
import logging
import os
import platform
import re
import subprocess
from collections import Counter
from pathlib import Path

import h5py
import yaml

from postprocessing.run_outputs import RunOutputs
from postprocessing.timeseries_model import seasonal_terms_supported

logger = logging.getLogger(__name__)

METADATA_KEY = "insar_processing"
SIDECAR_SUFFIX = ".metadata.json"
METADATA_FORMAT_VERSION = "1"

PIPELINE_CONFIG_NAME = "pipeline_config.yaml"
VERSIONS_FILE = Path("logs") / "versions.json"
ISCE_LOG_DIR = Path("logs") / "isce"
ISCE_UNPACK_LOG_PATTERN = "run_0[12]_*.log"
DOLPHIN_CONFIG = Path("dolphin") / "dolphin_config.yaml"
DOLPHIN_LOG = Path("dolphin") / "dolphin.log"
MINTPY_LOG = Path("logs") / "mintpy.log"
MINTPY_TEMPLATE_DEFAULT = "smallbaselineApp.cfg"

# dolphin logs JSON lines, so the version ends before a quote or comma.
DOLPHIN_VERSION_PATTERN = re.compile(r"dolphin version: (?P<version>[^\s\",]+)")
MINTPY_VERSION_PATTERN = re.compile(r"MintPy version (?P<version>[^,\s]+)")
DIRNAME_PATTERN = re.compile(r"'--dirname', '(?P<safes>[^']+)'")
ORBIT_FILE_PATTERN = re.compile(r"Extracting orbit from Orbit File:\s+(?P<orbit>\S+)")
SAFE_PATTERN = re.compile(
    r"(?P<mission>S1[A-D])_IW_SLC__\w{4}_(?P<date>\d{8})T\d{6}_\d{8}T\d{6}_\d{6}_\w{6}_\w{4}")
ORBIT_TYPES = ("POEORB", "RESORB")

RUN_ENVIRONMENT_PREFIX = "INSAR_RUN_"
GIT_COMMIT_VARIABLE = "POSTPROCESSING_GIT_COMMIT"
DIRECTION_ABBREVIATIONS = {"ascending": "ASC", "descending": "DES"}
FILE_NAME_SENSOR = "S1"
UNSAFE_NAME_CHARACTERS = re.compile(r"[^A-Za-z0-9-]+")
LIBRARIES = ("numpy", "pyarrow", "rasterio", "h5py", "pyproj", "shapely")


def _text(value):
    """Return an HDF5 attribute as str (None stays None)."""
    if value is None:
        return None
    return value.decode() if isinstance(value, bytes) else str(value)


def _read_yaml(path: Path):
    """Read a YAML file, or return None (with a warning) if missing or invalid."""
    if not path.is_file():
        logger.warning("%s not found; its settings are left out of the metadata.", path)
        return None
    try:
        return yaml.safe_load(path.read_text())
    except yaml.YAMLError as error:
        logger.warning("Cannot read %s: %s", path, error)
        return None


def read_mintpy_template(mintpy_dir: Path) -> dict | None:
    """Read the pipeline's MintPy template (``key = value`` lines) as a dict.

    The pipeline writes its own template next to MintPy's generated
    ``smallbaselineApp.cfg``; that one is used. Values stay strings.
    """
    templates = [path for path in sorted(mintpy_dir.glob("*.cfg"))
                 if path.name != MINTPY_TEMPLATE_DEFAULT]
    if not templates:
        return None
    options = {}
    for line in templates[0].read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            key, value = (part.strip() for part in line.split("=", 1))
            options[key] = value
    return options


def read_slc_orbits(isce_log_dir: Path) -> list[dict]:
    """Return each SLC with the orbit file ISCE2 used for it, from the unpack logs.

    Each ``run_01``/``run_02`` command log names the SAFE(s) (``--dirname``)
    followed by the orbit file ISCE2 extracted the orbit from.

    Returns
    -------
    list[dict]
        ``date``, ``mission``, ``slc``, ``orbit_files``, ``orbit_type`` per SAFE, by date.
    """
    slcs = {}
    for log_file in sorted(isce_log_dir.glob(ISCE_UNPACK_LOG_PATTERN)):
        current_safes = []
        for line in log_file.read_text(errors="replace").splitlines():
            dirname = DIRNAME_PATTERN.search(line)
            if dirname:
                current_safes = [Path(safe).name for safe in dirname["safes"].split(",")]
                for safe in current_safes:
                    slcs.setdefault(safe, {"orbit_files": []})
                continue
            orbit = ORBIT_FILE_PATTERN.search(line)
            if orbit:
                for safe in current_safes:
                    name = Path(orbit["orbit"]).name
                    if name not in slcs[safe]["orbit_files"]:
                        slcs[safe]["orbit_files"].append(name)
    result = []
    for safe, info in slcs.items():
        match = SAFE_PATTERN.search(safe)
        orbit_types = {kind for kind in ORBIT_TYPES for orbit in info["orbit_files"] if kind in orbit}
        result.append({
            "date": f"{match['date'][:4]}-{match['date'][4:6]}-{match['date'][6:]}" if match else None,
            "mission": match["mission"] if match else None,
            "slc": safe,
            "orbit_files": info["orbit_files"],
            "orbit_type": orbit_types.pop() if len(orbit_types) == 1 else None,
        })
    return sorted(result, key=lambda slc: (slc["date"] or "", slc["slc"]))


def _version_from_log(path: Path, pattern: re.Pattern) -> str | None:
    """Return the last version a log file reports, or None."""
    if not path.is_file():
        return None
    version = None
    with path.open(errors="replace") as lines:
        for line in lines:
            match = pattern.search(line)
            if match:
                version = match["version"]
    return version


def read_pipeline_versions(run_dir: Path) -> dict:
    """Return ISCE2, dolphin, MintPy and pipeline versions of the run.

    Uses ``logs/versions.json`` (written by the pipeline); for runs without it,
    dolphin and MintPy versions are taken from their logs.
    """
    versions_file = run_dir / VERSIONS_FILE
    if versions_file.is_file():
        return json.loads(versions_file.read_text())
    logger.warning("%s not found; taking dolphin and MintPy versions from their logs.",
                   versions_file)
    return {
        "isce2": None,
        "dolphin": _version_from_log(run_dir / DOLPHIN_LOG, DOLPHIN_VERSION_PATTERN),
        "mintpy": _version_from_log(run_dir / MINTPY_LOG, MINTPY_VERSION_PATTERN),
        "isce_dolphin_mintpy": None,
    }


def postprocessing_version() -> str | None:
    """Return the git commit of this code (``-dirty`` with uncommitted changes), or None."""
    if os.environ.get(GIT_COMMIT_VARIABLE):
        return os.environ[GIT_COMMIT_VARIABLE]
    repository = Path(__file__).resolve().parent.parent
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repository,
                                capture_output=True, text=True, check=True).stdout.strip()
        changes = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"],
                                 cwd=repository, capture_output=True, text=True,
                                 check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None
    return f"{commit}-dirty" if changes else commit


def _library_versions() -> dict:
    """Return the versions of Python and the libraries used for the export."""
    versions = {"python": platform.python_version()}
    for library in LIBRARIES:
        try:
            versions[library] = importlib.metadata.version(library)
        except importlib.metadata.PackageNotFoundError:
            versions[library] = None
    return versions


def _reference_point(outputs: RunOutputs, attributes: dict) -> dict:
    """Return MintPy's reference pixel with its latitude and longitude."""
    ref_y, ref_x = attributes.get("REF_Y"), attributes.get("REF_X")
    if ref_y is None or ref_x is None:
        return {"row": None, "col": None, "lat": None, "lon": None}
    row, col = int(ref_y), int(ref_x)
    position = {"lat": None, "lon": None}
    for name in ("lat", "lon"):
        layer = outputs.layers.get(name)
        if layer is not None:
            with h5py.File(layer.path, "r") as h5_file:
                position[name] = float(h5_file[layer.dataset][row, col])
    return {"row": row, "col": col, **position}


def _swaths(run_dir: Path, attributes: dict) -> list[str]:
    """Return the processed swaths, e.g. ["IW1", "IW2"]."""
    reference_dir = run_dir / "isce" / "reference"
    swaths = sorted(path.stem for path in reference_dir.glob("IW[1-3].xml"))
    if swaths:
        return swaths
    numbers = attributes.get("swathNumber") or ""
    return [f"IW{number}" for number in str(numbers) if number.isdigit()]


def collect_metadata(outputs: RunOutputs, area_name: str | None = None,
                     wavelength_m: float | None = None, excluded: list[str] | None = None,
                     all_pixels: bool = False) -> dict:
    """Collect everything known about how this run was made.

    Returns
    -------
    dict
        JSON-serializable metadata with the sections ``product``, ``software``,
        ``inputs``, ``settings`` and ``run``. Unknown values are None.
    """
    run_dir = outputs.run_dir
    with h5py.File(outputs.timeseries_file, "r") as h5_file:
        attributes = {key: _text(value) for key, value in h5_file.attrs.items()}
    pipeline_config = _read_yaml(run_dir / PIPELINE_CONFIG_NAME) or {}
    paths = pipeline_config.get("paths") or {}
    slcs = read_slc_orbits(run_dir / ISCE_LOG_DIR)
    used_dates = {date.isoformat() for date in outputs.dates}
    slcs = [slc for slc in slcs if slc["date"] in used_dates or slc["date"] is None]
    track = attributes.get("trackNumber")
    dates = outputs.dates

    return {
        "metadata_format_version": METADATA_FORMAT_VERSION,
        "product": {
            "sensor": "Sentinel-1",
            "missions": sorted({slc["mission"] for slc in slcs if slc["mission"]}),
            "track": int(track) if track not in (None, "") else None,
            "orbit_direction": attributes.get("ORBIT_DIRECTION"),
            "swaths": _swaths(run_dir, attributes),
            "polarization": (attributes.get("POLARIZATION") or "").upper() or None,
            "coordinates": "radar grid; one point per pixel at its ISCE2 latitude/longitude",
            "area_name": area_name,
            "first_date": dates[0].isoformat(),
            "last_date": dates[-1].isoformat(),
            "reference_date": outputs.reference_date.isoformat(),
            "reference_point": _reference_point(outputs, attributes),
            "number_of_dates": len(dates),
            "dates": [date.isoformat() for date in dates],
        },
        "software": {
            **read_pipeline_versions(run_dir),
            "postprocessing": postprocessing_version(),
            **{f"postprocessing_{key}": value for key, value in _library_versions().items()},
        },
        "inputs": {
            "slcs": slcs,
            "orbit_types": dict(Counter(slc["orbit_type"] for slc in slcs if slc["orbit_type"])),
            "dem_file": paths.get("dem_file"),
            "dem_vertical_datum": paths.get("dem_vertical_datum"),
            "water_mask_file": (pipeline_config.get("water_mask") or {}).get("file"),
        },
        "settings": {
            "pipeline": pipeline_config or None,
            "dolphin": _read_yaml(run_dir / DOLPHIN_CONFIG),
            "mintpy_template": read_mintpy_template(outputs.mintpy_dir),
            "mintpy_timeseries_file": outputs.timeseries_file.name,
            "radar_grid": {
                "rows": outputs.shape[0], "columns": outputs.shape[1],
                "range_pixel_size_m": _float(attributes.get("RANGE_PIXEL_SIZE")),
                "azimuth_pixel_size_m": _float(attributes.get("AZIMUTH_PIXEL_SIZE")),
            },
            "postprocessing": {
                "wavelength_m": wavelength_m,
                "excluded": sorted(excluded or []),
                "all_pixels": all_pixels,
                "seasonal_terms_fitted": bool(seasonal_terms_supported(dates)),
            },
        },
        "run": {
            "processed_at_utc": datetime.datetime.now(datetime.timezone.utc)
            .isoformat(timespec="seconds"),
            **{name.removeprefix(RUN_ENVIRONMENT_PREFIX).lower(): value
               for name, value in sorted(os.environ.items())
               if name.startswith(RUN_ENVIRONMENT_PREFIX)},
        },
    }


def _float(value):
    """Convert a numeric attribute to float, or None."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def clean_area_name(area_name: str) -> str:
    """Make an area name safe for file names, e.g. ``"Wien Nord"`` -> ``"Wien-Nord"``.

    Raises
    ------
    ValueError
        If nothing usable is left.
    """
    cleaned = UNSAFE_NAME_CHARACTERS.sub("-", area_name).strip("-")
    if not cleaned:
        raise ValueError(f"Area name {area_name!r} contains no letters or digits.")
    return cleaned


def build_output_name(metadata: dict) -> str:
    """Build the file name from the metadata.

    Pattern: ``S1_<ASC|DES>_<track>_<swaths>_<pol>_<first>_<last>_<area>.parquet``,
    e.g. ``S1_ASC_073_IW1_VV_20220101_20220125_vienna.parquet`` (same as
    dolphin_postprocessing). Unknown parts are left out with a warning.
    """
    product = metadata["product"]
    direction = product["orbit_direction"]
    track = product["track"]
    parts = {
        "orbit direction": DIRECTION_ABBREVIATIONS.get(direction.lower())
        if isinstance(direction, str) else None,
        "track": f"{track:03d}" if isinstance(track, int) else None,
        "subswaths": "-".join(product["swaths"]) or None,
        "polarization": product["polarization"],
    }
    missing = [description for description, part in parts.items() if part is None]
    if missing:
        logger.warning("Unknown %s: left out of the file name.", ", ".join(missing))
    name_parts = [FILE_NAME_SENSOR, *(part for part in parts.values() if part)]
    name_parts += [product["first_date"].replace("-", ""), product["last_date"].replace("-", "")]
    if product["area_name"]:
        name_parts.append(clean_area_name(product["area_name"]))
    return "_".join(name_parts) + ".parquet"


def write_sidecar(metadata: dict, output_file: Path) -> Path:
    """Write the metadata as readable JSON next to the point file; return its path."""
    sidecar = output_file.with_name(output_file.stem + SIDECAR_SUFFIX)
    sidecar.write_text(json.dumps(metadata, indent=2, default=str) + "\n")
    return sidecar
