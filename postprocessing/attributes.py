"""The columns of the point file: names, groups, units and descriptions.

Every column belongs to a group. ``--exclude`` takes group names or single column
names, so for example ``--exclude displacement_std`` leaves out all per-date
standard deviations. ``point_id`` and the point geometry are always written.

The per-date columns are generated from the acquisition dates:
``D<YYYYMMDD>`` (displacement) and ``D<YYYYMMDD>_std`` (its standard deviation).
"""

import datetime
from dataclasses import dataclass

import pyarrow as pa


@dataclass(frozen=True)
class Attribute:
    """One column of the point file.

    Attributes
    ----------
    name : str
        Column name.
    group : str
        Group name, usable with ``--exclude``.
    data_type : pyarrow.DataType
        Column type in the Parquet file.
    unit : str
        Unit, or an empty string for unitless values.
    description : str
        What the value means, stored in the Parquet file.
    """

    name: str
    group: str
    data_type: pa.DataType
    unit: str
    description: str


# Group names.
POSITION = "position"
GEOMETRY = "geometry"
WATER = "water"
MINTPY_QUALITY = "mintpy_quality"
PHASE_LINKING = "phase_linking"
PERSISTENT_SCATTERERS = "persistent_scatterers"
SUBPIXEL = "subpixel"
MINTPY_VELOCITY = "mintpy_velocity"
DEM_ERROR = "dem_error"
MODEL = "model"
MODEL_COHERENCE = "model_coherence"
DISPLACEMENT = "displacement"
DISPLACEMENT_STD = "displacement_std"

POINT_ID = "point_id"

FIXED_ATTRIBUTES = (
    Attribute(POINT_ID, POSITION, pa.int64(), "",
              "Raster row * raster width + raster column. Stable for one output grid."),
    Attribute("row", POSITION, pa.int32(), "pixel", "Row in the radar grid (azimuth)."),
    Attribute("col", POSITION, pa.int32(), "pixel", "Column in the radar grid (range)."),
    Attribute("lon", POSITION, pa.float64(), "deg",
              "Longitude (WGS84) of the pixel center, from the ISCE2 geometry."),
    Attribute("lat", POSITION, pa.float64(), "deg",
              "Latitude (WGS84) of the pixel center, from the ISCE2 geometry."),
    Attribute("height", GEOMETRY, pa.float32(), "m",
              "DEM height above the WGS84 ellipsoid, sampled at the pixel center."),
    Attribute("incidence_angle", GEOMETRY, pa.float32(), "deg",
              "Angle between the line of sight and the ellipsoid normal."),
    Attribute("azimuth_angle", GEOMETRY, pa.float32(), "deg",
              "Azimuth of the ground-to-satellite line of sight, from north, anti-clockwise."),
    Attribute("los_east", GEOMETRY, pa.float32(), "",
              "East component of the unit vector from ground to satellite."),
    Attribute("los_north", GEOMETRY, pa.float32(), "",
              "North component of the unit vector from ground to satellite."),
    Attribute("los_up", GEOMETRY, pa.float32(), "",
              "Up component of the unit vector from ground to satellite."),
    Attribute("shadow_layover", GEOMETRY, pa.int8(), "",
              "1 if ISCE2 marks the pixel as shadow or layover, else 0."),
    Attribute("water_fraction", WATER, pa.float32(), "",
              "Fraction of the pixel area that is water in the water mask (0..1)."),
    Attribute("is_water", WATER, pa.int8(), "",
              "1 if the pipeline masked the pixel as water (fraction above its threshold)."),
    Attribute("temporal_coherence", MINTPY_QUALITY, pa.float32(), "",
              "MintPy temporal coherence of the network inversion (dolphin's phase "
              "linking coherence when dolphin inverted the network)."),
    Attribute("avg_spatial_coherence", MINTPY_QUALITY, pa.float32(), "",
              "Average interferometric coherence over all interferograms."),
    Attribute("mask_temporal_coherence", MINTPY_QUALITY, pa.int8(), "",
              "1 if the pixel is in MintPy's final mask (maskTempCoh.h5)."),
    Attribute("conncomp_ok", MINTPY_QUALITY, pa.int8(), "",
              "1 if the pixel is in a valid connected component (maskConnComp.h5)."),
    Attribute("num_inverted_ifgrams", MINTPY_QUALITY, pa.float32(), "",
              "Number of interferograms used for this pixel in the network inversion."),
    Attribute("phase_linking_coherence", PHASE_LINKING, pa.float32(), "",
              "dolphin's phase linking temporal coherence."),
    Attribute("similarity", PHASE_LINKING, pa.float32(), "",
              "dolphin's phase similarity to neighboring pixels."),
    Attribute("shp_count", PHASE_LINKING, pa.int32(), "",
              "Statistically homogeneous pixels used in dolphin's estimation window."),
    Attribute("is_ps", PERSISTENT_SCATTERERS, pa.int8(), "",
              "1 if at least one full-resolution PS pixel lies in this output pixel."),
    Attribute("amp_dispersion", PERSISTENT_SCATTERERS, pa.float32(), "",
              "Lowest amplitude dispersion in this output pixel. Lower is more PS-like."),
    # Optional: only written when the pipeline's subpixel step ran (subpixel/ exists).
    Attribute("lon_subpixel", SUBPIXEL, pa.float64(), "deg",
              "Longitude (WGS84) of the PS at its sub-pixel position, on the DEM surface "
              "(dem_error not applied). Empty where there is no PS."),
    Attribute("lat_subpixel", SUBPIXEL, pa.float64(), "deg",
              "Latitude (WGS84) of the PS at its sub-pixel position, on the DEM surface "
              "(dem_error not applied). Empty where there is no PS."),
    Attribute("height_subpixel", SUBPIXEL, pa.float32(), "m",
              "DEM height above the WGS84 ellipsoid at the sub-pixel position."),
    Attribute("shift_east", SUBPIXEL, pa.float32(), "m",
              "Eastward shift from the pixel center (lon, lat) to the sub-pixel position. "
              "Exact with float64 lat/lon (current pipeline); older runs stored float32, "
              "~0.2 m rounding."),
    Attribute("shift_north", SUBPIXEL, pa.float32(), "m",
              "Northward shift from the pixel center (lon, lat) to the sub-pixel position. "
              "Exact with float64 lat/lon (current pipeline); older runs stored float32, "
              "~0.2 m rounding."),
    Attribute("dem_error_subpixel", SUBPIXEL, pa.float32(), "m",
              "dem_error re-referenced to the sub-pixel position: dem_error - range_offset * "
              "full-resolution slant-range spacing * cos(incidence_angle). Add this, not dem_error, "
              "to height_subpixel; adding dem_error would count the range offset twice."),
    Attribute("range_offset", SUBPIXEL, pa.float32(), "pixel",
              "Sub-pixel offset in full-resolution range samples. Positive = farther from the sensor."),
    Attribute("azimuth_offset", SUBPIXEL, pa.float32(), "pixel",
              "Sub-pixel offset in full-resolution azimuth lines. Positive = later azimuth time."),
    Attribute("range_offset_std", SUBPIXEL, pa.float32(), "pixel",
              "Standard deviation of the per-date range offsets. Large = several scatterers or "
              "noise; do not use the offset then."),
    Attribute("azimuth_offset_std", SUBPIXEL, pa.float32(), "pixel",
              "Standard deviation of the per-date azimuth offsets. Large = several scatterers or "
              "noise; do not use the offset then."),
    Attribute("subpixel_valid_dates", SUBPIXEL, pa.int16(), "",
              "Number of dates used for the sub-pixel estimate."),
    Attribute("subpixel_on_peak", SUBPIXEL, pa.int8(), "",
              "1 = this PS pixel is the amplitude peak itself; 0 = it lies on the flank of a peak "
              "in a neighbouring pixel (same scatterer). Keep 1 for one point per scatterer."),
    Attribute("velocity_mintpy", MINTPY_VELOCITY, pa.float32(), "mm/year",
              "MintPy's velocity: linear fit of the final, corrected time series."),
    Attribute("velocity_mintpy_std", MINTPY_VELOCITY, pa.float32(), "mm/year",
              "Standard deviation of velocity_mintpy."),
    Attribute("dem_error", DEM_ERROR, pa.float32(), "m",
              "DEM error estimated by MintPy's topographic residual correction."),
    Attribute("mean_velocity", MODEL, pa.float32(), "mm/year",
              "EGMS: linear term of a linear + seasonal model."),
    Attribute("mean_velocity_std", MODEL, pa.float32(), "mm/year",
              "EGMS: standard deviation of mean_velocity."),
    Attribute("acceleration", MODEL, pa.float32(), "mm/year^2",
              "EGMS: twice the quadratic coefficient of a quadratic + seasonal model."),
    Attribute("acceleration_std", MODEL, pa.float32(), "mm/year^2",
              "EGMS: standard deviation of acceleration."),
    Attribute("seasonality", MODEL, pa.float32(), "mm",
              "EGMS: amplitude of the annual term of a cubic + seasonal model."),
    Attribute("seasonality_std", MODEL, pa.float32(), "mm",
              "EGMS: standard deviation of seasonality (Rayleigh distribution)."),
    Attribute("season_peak_doy", MODEL, pa.float32(), "day",
              "Day of the year (0 = 1 January) when the annual term peaks toward the satellite."),
    Attribute("season_peak_doy_std", MODEL, pa.float32(), "day",
              "Standard deviation of season_peak_doy. Unreliable when seasonality is small."),
    Attribute("rmse", MODEL, pa.float32(), "mm",
              "EGMS: root mean square of the cubic + seasonal model residuals."),
    Attribute("model_temporal_coherence", MODEL_COHERENCE, pa.float32(), "",
              "EGMS-style temporal coherence of the linear + seasonal model residuals."),
)

GROUPS = tuple(dict.fromkeys(
    [attribute.group for attribute in FIXED_ATTRIBUTES]
    + [DISPLACEMENT, DISPLACEMENT_STD]
))
ALWAYS_WRITTEN = {POINT_ID}
# Sub-pixel shift in meters, computed from the pixel center and the sub-pixel position.
SHIFT_COLUMNS = ("shift_east", "shift_north")
# DEM error at the sub-pixel position, computed from dem_error, range_offset and incidence_angle.
DEM_ERROR_SUBPIXEL = "dem_error_subpixel"
# Columns that need the time-series model fit.
MODEL_COLUMNS = frozenset(
    attribute.name
    for attribute in FIXED_ATTRIBUTES
    if attribute.group in (MODEL, MODEL_COHERENCE)
)


def displacement_column(date: datetime.date) -> str:
    """Return the displacement column name for a date, e.g. ``D20230108``."""
    return f"D{date:%Y%m%d}"


def displacement_std_column(date: datetime.date) -> str:
    """Return the standard deviation column name for a date."""
    return f"{displacement_column(date)}_std"


def date_attributes(dates: list[datetime.date]) -> list[Attribute]:
    """Create the per-date columns.

    Parameters
    ----------
    dates : list[datetime.date]
        All acquisition dates of the MintPy time series.

    Returns
    -------
    list[Attribute]
        All displacement columns, then all standard deviation columns.
    """
    displacement = [
        Attribute(displacement_column(date), DISPLACEMENT, pa.float32(), "mm",
                  "Line-of-sight displacement since the reference date, after MintPy's "
                  "corrections, relative to the reference point. Positive toward the satellite.")
        for date in dates
    ]
    std = [
        Attribute(displacement_std_column(date), DISPLACEMENT_STD, pa.float32(), "mm",
                  "Displacement standard deviation from dolphin's CRLB. Phase linking "
                  "only: no unwrapping, atmospheric or correction errors.")
        for date in dates
    ]
    return displacement + std


def select_attributes(
    dates: list[datetime.date], excluded: list[str] | None = None
) -> list[Attribute]:
    """Return the columns to write, without the excluded groups and columns.

    Parameters
    ----------
    dates : list[datetime.date]
        All acquisition dates.
    excluded : list[str], optional
        Group names and fixed column names to leave out.

    Returns
    -------
    list[Attribute]
        The selected columns, in output order.

    Raises
    ------
    ValueError
        If a name is neither a group nor a fixed column, or cannot be excluded.
    """
    excluded = set(excluded or [])
    fixed_names = {attribute.name for attribute in FIXED_ATTRIBUTES}

    unknown = sorted(excluded - set(GROUPS) - fixed_names)
    if unknown:
        raise ValueError(
            f"Cannot exclude {unknown}: not a group or column. "
            f"Groups: {list(GROUPS)}. Columns: {sorted(fixed_names)}."
        )
    always_written = sorted(excluded & ALWAYS_WRITTEN)
    if always_written:
        raise ValueError(f"{always_written} are always written.")

    return [
        attribute
        for attribute in [*FIXED_ATTRIBUTES, *date_attributes(dates)]
        if attribute.name not in excluded and attribute.group not in excluded
    ]
