"""Write the results of an isce2-dolphin-mintpy run as one GeoParquet point file.

Every selected pixel of the radar grid becomes one point at its own latitude and
longitude (from the ISCE2 geometry). By default only pixels in MintPy's final
mask (``maskTempCoh.h5``) are written; ``all_pixels`` writes every pixel with a
valid position, with the mask as a column to filter later.

The grid is read in blocks of rows, so the whole area never has to fit in
memory. The EGMS-style model is fitted in the same loop. The output goes to
``<output>.partial`` first and is renamed at the end, so an interrupted run
never leaves a file that looks complete.
"""

import json
import logging
import math
import warnings
from contextlib import ExitStack
from pathlib import Path

import h5py
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import rasterio
import shapely
from pyproj import CRS
from rasterio.errors import NotGeoreferencedWarning
from rasterio.windows import Window

from postprocessing.attributes import (
    DISPLACEMENT,
    DISPLACEMENT_STD,
    MODEL,
    MODEL_COHERENCE,
    MODEL_COLUMNS,
    SHIFT_COLUMNS,
    SUBPIXEL,
    Attribute,
    displacement_column,
    displacement_std_column,
    select_attributes,
)
from postprocessing.run_outputs import Layer, RunOutputs
from postprocessing.timeseries_model import estimate_model_parameters, seasonal_terms_supported

logger = logging.getLogger(__name__)

DEFAULT_BLOCK_ROWS = 256
GEOMETRY_COLUMN = "geometry"
GEOPARQUET_VERSION = "1.1.0"
PARQUET_COMPRESSION = "zstd"
# Points are stored as (longitude, latitude); OGC:CRS84 states that axis order.
POINT_CRS = "OGC:CRS84"
MM_PER_M = 1000.0

# Columns read from one layer, converted with a factor (m -> mm etc.).
LAYER_SCALE = {"velocity_mintpy": MM_PER_M, "velocity_mintpy_std": MM_PER_M}
LOS_SOURCES = ("incidence_angle", "azimuth_angle")
LOS_COLUMNS = ("los_east", "los_north", "los_up")
# Not written as columns themselves (lat/lon become position columns).
POSITION_LAYERS = ("lat", "lon")
# Layers the sub-pixel shift is computed from (besides lat/lon).
SHIFT_SOURCES = ("lat_subpixel", "lon_subpixel")

# WGS84 ellipsoid, for converting small latitude/longitude differences to meters.
WGS84_SEMI_MAJOR_AXIS_M = 6378137.0
WGS84_ECCENTRICITY_SQUARED = 6.69437999014e-3


class ExportError(RuntimeError):
    """The point file cannot be written. The message explains why."""


class _LayerReader:
    """Reads rows of HDF5 datasets and GeoTIFFs, keeping the files open."""

    def __init__(self, stack: ExitStack):
        self._stack = stack
        self._open = {}
        # dolphin's rasters are in radar coordinates and have no geotransform on
        # purpose; rasterio would warn about that for every file.
        warnings.filterwarnings("ignore", category=NotGeoreferencedWarning)

    def _handle(self, path: Path, is_hdf5: bool):
        if path not in self._open:
            opener = h5py.File(path, "r") if is_hdf5 else rasterio.open(path)
            self._open[path] = self._stack.enter_context(opener)
        return self._open[path]

    def shape(self, layer: Layer) -> tuple[int, int]:
        """Return the (rows, columns) of a layer."""
        if layer.dataset is not None:
            return tuple(self._handle(layer.path, True)[layer.dataset].shape[-2:])
        return self._handle(layer.path, False).shape

    def rows(self, layer: Layer, first: int, last: int) -> np.ndarray:
        """Read rows ``first..last-1`` as float64, NaN for no-data."""
        if layer.dataset is not None:
            return self._handle(layer.path, True)[layer.dataset][first:last].astype(np.float64)
        source = self._handle(layer.path, False)
        window = Window(0, first, source.width, last - first)
        return source.read(1, window=window, masked=True).astype(np.float64).filled(np.nan)

    def timeseries_rows(self, path: Path, first: int, last: int) -> np.ndarray:
        """Read rows of all dates of a MintPy time series: shape (dates, rows, cols)."""
        return self._handle(path, True)["timeseries"][:, first:last].astype(np.float64)


def _usable_layers(outputs: RunOutputs, reader: _LayerReader) -> dict[str, Layer]:
    """Return the layers on the time series grid; warn about the others."""
    usable = {}
    for name, layer in outputs.layers.items():
        shape = reader.shape(layer)
        if shape == outputs.shape:
            usable[name] = layer
        else:
            logger.warning("Leaving out %s: %s has shape %s, the time series %s.",
                           name, layer.path.name, shape, outputs.shape)
    return usable


def _crlb_usable(outputs: RunOutputs, reader: _LayerReader) -> bool:
    """True if every date has a CRLB raster on the time series grid."""
    return bool(outputs.crlb_files) and all(
        date in outputs.crlb_files and reader.shape(Layer(outputs.crlb_files[date])) == outputs.shape
        for date in outputs.dates
    )


def _without_unused_subpixel(selected: list[Attribute], layers: dict[str, Layer]) -> list[Attribute]:
    """Drop the sub-pixel columns if the run has no sub-pixel output at all.

    The subpixel step is optional and off by default, so its absence is normal
    and not worth a warning. If only some of its files are missing, the usual
    warning in ``_available_attributes`` names them.
    """
    subpixel_names = {attribute.name for attribute in selected if attribute.group == SUBPIXEL}
    if not subpixel_names or subpixel_names & set(layers):
        return selected
    logger.info("No sub-pixel output (subpixel step not run); sub-pixel columns not written.")
    return [attribute for attribute in selected if attribute.group != SUBPIXEL]


def _available_attributes(selected: list[Attribute], layers: dict[str, Layer],
                          crlb_usable: bool) -> list[Attribute]:
    """Drop the selected columns whose sources are missing; warn about them."""
    selected = _without_unused_subpixel(selected, layers)

    def is_available(attribute: Attribute) -> bool:
        if attribute.name in LOS_COLUMNS:
            return all(source in layers for source in LOS_SOURCES)
        if attribute.name in SHIFT_COLUMNS:
            return all(source in layers for source in SHIFT_SOURCES)
        if attribute.group == DISPLACEMENT_STD:
            return crlb_usable
        if attribute.group in (DISPLACEMENT, MODEL, MODEL_COHERENCE) or attribute.name in (
                "point_id", "row", "col", "lon", "lat"):
            return True
        return attribute.name in layers

    available = [attribute for attribute in selected if is_available(attribute)]
    dropped = [attribute.name for attribute in selected if attribute not in available]
    if dropped:
        logger.warning("Sources missing or unusable. Leaving out: %s", ", ".join(dropped))
    return available


def _build_schema(attributes: list[Attribute], file_metadata: dict[str, str] | None) -> pa.Schema:
    """Build the Parquet schema with units, descriptions and GeoParquet metadata."""
    fields = [
        pa.field(attribute.name, attribute.data_type,
                 metadata={"unit": attribute.unit, "description": attribute.description})
        for attribute in attributes
    ]
    fields.append(pa.field(GEOMETRY_COLUMN, pa.binary()))
    geo_metadata = {
        "version": GEOPARQUET_VERSION,
        "primary_column": GEOMETRY_COLUMN,
        "columns": {
            GEOMETRY_COLUMN: {
                "encoding": "WKB",
                "geometry_types": ["Point"],
                "crs": CRS.from_string(POINT_CRS).to_json_dict(),
            }
        },
    }
    return pa.schema(fields, metadata={**(file_metadata or {}), "geo": json.dumps(geo_metadata)})


def _to_arrow_array(values: np.ndarray, data_type: pa.DataType) -> pa.Array:
    """Convert values with NaN for missing data to an Arrow array with nulls."""
    missing = np.isnan(values)
    if pa.types.is_integer(data_type):
        values = np.where(missing, 0, values)
    return pa.array(values.astype(data_type.to_pandas_dtype()), type=data_type, mask=missing)


def los_unit_vector(incidence_deg: np.ndarray, azimuth_deg: np.ndarray):
    """Return (east, north, up) of the unit vector from ground to satellite.

    MintPy/ISCE2 convention: ``azimuth_deg`` is the azimuth of that vector,
    measured from north, anti-clockwise (same as MintPy's
    ``get_unit_vector4component_of_interest``).
    """
    incidence = np.deg2rad(incidence_deg)
    azimuth = np.deg2rad(azimuth_deg)
    return (-np.sin(incidence) * np.sin(azimuth),
            np.sin(incidence) * np.cos(azimuth),
            np.cos(incidence))


def shift_in_meters(lat_deg: np.ndarray, lon_deg: np.ndarray,
                    target_lat_deg: np.ndarray, target_lon_deg: np.ndarray):
    """Return the (east, north) distance in meters from one position to a nearby one.

    Uses the WGS84 radii of curvature at ``lat_deg`` (meridian radius for north,
    prime vertical radius for east). For shifts of a few meters the error of
    this local approximation is far below a millimeter.

    Parameters
    ----------
    lat_deg, lon_deg : np.ndarray
        Start positions (pixel centers), degrees.
    target_lat_deg, target_lon_deg : np.ndarray
        End positions (sub-pixel positions), degrees. NaN gives NaN.

    Returns
    -------
    tuple of np.ndarray
        Eastward and northward shift in meters.
    """
    lat = np.deg2rad(lat_deg)
    denominator = 1.0 - WGS84_ECCENTRICITY_SQUARED * np.sin(lat) ** 2
    meridian_radius = (WGS84_SEMI_MAJOR_AXIS_M * (1.0 - WGS84_ECCENTRICITY_SQUARED)
                       / denominator ** 1.5)
    prime_vertical_radius = WGS84_SEMI_MAJOR_AXIS_M / np.sqrt(denominator)
    east = np.deg2rad(target_lon_deg - lon_deg) * prime_vertical_radius * np.cos(lat)
    north = np.deg2rad(target_lat_deg - lat_deg) * meridian_radius
    return east, north


def _selected_pixels(block: dict[str, np.ndarray], all_pixels: bool) -> np.ndarray:
    """Boolean mask of the pixels to write in one block of rows.

    A pixel needs a position (ISCE2 marks no-data with lat = lon = 0). Unless
    ``all_pixels`` is set, it must also be in MintPy's final mask.
    """
    lat, lon = block["lat"], block["lon"]
    selected = np.isfinite(lat) & np.isfinite(lon) & ~((lat == 0) & (lon == 0))
    if not all_pixels:
        selected &= block["mask_temporal_coherence"] > 0
    return selected


def export_points(
    outputs: RunOutputs,
    output_file: Path,
    wavelength_m: float,
    excluded: list[str] | None = None,
    all_pixels: bool = False,
    block_rows: int = DEFAULT_BLOCK_ROWS,
    overwrite: bool = False,
    file_metadata: dict[str, str] | None = None,
) -> int:
    """Write the selected pixels of a run as points to GeoParquet.

    Parameters
    ----------
    outputs : RunOutputs
        From ``find_run_outputs``.
    output_file : Path
        The ``.parquet`` file to write.
    wavelength_m : float
        Radar wavelength in meters (CRLB and model coherence columns).
    excluded : list[str], optional
        Groups and columns to leave out. See ``attributes.py``.
    all_pixels : bool, optional
        Write every pixel with a position, not only those in MintPy's mask.
    block_rows : int, optional
        Radar rows read at a time.
    overwrite : bool, optional
        Replace ``output_file`` if it exists.
    file_metadata : dict[str, str], optional
        Extra file-level metadata entries (key -> text).

    Returns
    -------
    int
        Number of points written.

    Raises
    ------
    ExportError
        If the output exists and ``overwrite`` is False, MintPy's mask is
        missing without ``all_pixels``, or no pixel is selected.
    ValueError
        If ``excluded`` contains an unknown name.
    """
    if output_file.exists() and not overwrite:
        raise ExportError(f"{output_file} exists. Use --overwrite to replace it.")

    partial_file = output_file.with_name(f"{output_file.name}.partial")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    point_count = 0
    try:
        with ExitStack() as stack:
            reader = _LayerReader(stack)
            layers = _usable_layers(outputs, reader)
            if not all_pixels and "mask_temporal_coherence" not in layers:
                raise ExportError(
                    "MintPy's mask (maskTempCoh.h5) is missing. Use --all-pixels to "
                    "write every pixel without it.")
            attributes = _available_attributes(
                select_attributes(outputs.dates, excluded), layers, _crlb_usable(outputs, reader))
            names = {attribute.name for attribute in attributes}
            seasonal = seasonal_terms_supported(outputs.dates)
            if names & MODEL_COLUMNS and not seasonal:
                logger.warning("The dates span less than one year. Fitting the models without "
                               "seasonal terms; the seasonal columns stay empty.")

            schema = _build_schema(attributes, file_metadata)
            writer = stack.enter_context(
                pq.ParquetWriter(partial_file, schema, compression=PARQUET_COMPRESSION))
            needed_layers = {name for name in layers if name in names or name in POSITION_LAYERS
                             or name in LOS_SOURCES or name == "mask_temporal_coherence"
                             or (name in SHIFT_SOURCES and names & set(SHIFT_COLUMNS))}
            radians_to_mm = wavelength_m / (4 * math.pi) * MM_PER_M
            length, width = outputs.shape

            for first in range(0, length, block_rows):
                last = min(first + block_rows, length)
                block = {name: reader.rows(layers[name], first, last) for name in needed_layers}
                selected = _selected_pixels(block, all_pixels)
                if not selected.any():
                    continue
                rows, cols = np.nonzero(selected)
                columns = {name: values[rows, cols] for name, values in block.items()}
                for name, factor in LAYER_SCALE.items():
                    if name in columns:
                        columns[name] = columns[name] * factor
                if "is_water" in columns:
                    # waterMask is True for land; the column is 1 for water.
                    columns["is_water"] = 1.0 - columns["is_water"]
                if all(source in columns for source in LOS_SOURCES):
                    east, north, up = los_unit_vector(columns["incidence_angle"],
                                                      columns["azimuth_angle"])
                    columns |= {"los_east": east, "los_north": north, "los_up": up}
                if all(source in columns for source in SHIFT_SOURCES):
                    shift_east, shift_north = shift_in_meters(
                        columns["lat"], columns["lon"],
                        columns["lat_subpixel"], columns["lon_subpixel"])
                    columns |= {"shift_east": shift_east, "shift_north": shift_north}

                displacement_mm = (reader.timeseries_rows(outputs.timeseries_file, first, last)
                                   [:, rows, cols].T * MM_PER_M)
                for date, values in zip(outputs.dates, displacement_mm.T):
                    columns[displacement_column(date)] = values
                if names & MODEL_COLUMNS:
                    columns |= estimate_model_parameters(outputs.dates, displacement_mm,
                                                         wavelength_m, seasonal)
                if names & {displacement_std_column(date) for date in outputs.dates}:
                    for date in outputs.dates:
                        crlb = reader.rows(Layer(outputs.crlb_files[date]), first, last)
                        columns[displacement_std_column(date)] = crlb[rows, cols] * radians_to_mm

                absolute_rows = rows + first
                columns |= {"point_id": absolute_rows.astype(np.int64) * width + cols,
                            "row": absolute_rows, "col": cols}
                arrays = [_to_arrow_array(np.asarray(columns[a.name], dtype=np.float64), a.data_type)
                          for a in attributes]
                points = shapely.points(columns["lon"], columns["lat"])
                arrays.append(pa.array(shapely.to_wkb(points), pa.binary()))
                writer.write_table(pa.Table.from_arrays(arrays, schema=schema))
                point_count += rows.size
                logger.info("Rows %d-%d of %d: %d points so far.", first, last - 1, length,
                            point_count)

        if point_count == 0:
            raise ExportError("No pixel selected. Check MintPy's mask or use --all-pixels.")
        partial_file.replace(output_file)
    finally:
        partial_file.unlink(missing_ok=True)

    logger.info("Wrote %d points with %d columns to %s", point_count, len(attributes) + 1,
                output_file)
    return point_count
