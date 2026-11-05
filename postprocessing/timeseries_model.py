"""Fit EGMS-style models to displacement time series.

The definitions follow the EGMS Product Description and Format Specification
(sections 5.1 and 11.4). EGMS fits three separate least-squares models, with time
``t`` in years since the first acquisition:

- ``mean_velocity``: linear + seasonal model, value = linear coefficient.
- ``acceleration``: quadratic + seasonal model, value = twice the quadratic
  coefficient. Polynomial columns are written as ``t**k / k!`` here, so the
  coefficient of ``t**2 / 2`` is the acceleration itself.
- ``seasonality`` and ``rmse``: cubic + seasonal model. Seasonality is the
  amplitude of the annual sine and cosine terms.

EGMS has no seasonal phase. ``season_peak_doy`` is added here: the seasonal terms
use calendar time, so the phase converts to the day of the year at which the
annual term reaches its maximum (toward the satellite).

All series of one block are fitted at once. Pixels without gaps share one design
matrix, so they are solved with a single matrix product. Pixels with gaps are
solved one normal-equation system per pixel, vectorized with numpy.
"""

import datetime
import logging
import math
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger(__name__)

DAYS_PER_YEAR = 365.25
# The seasonal terms use time since this date, so their phase means a calendar day.
SEASONAL_EPOCH = datetime.date(2000, 1, 1)
# With less than one year of data, the seasonal terms cannot be separated from the
# polynomial, so they are left out of all models.
MINIMUM_SEASONAL_SPAN_YEARS = 1.0
# Standard deviation of a Rayleigh-distributed amplitude, relative to the
# standard deviation of its two Gaussian components (EGMS section 11.4).
RAYLEIGH_STD_FACTOR = math.sqrt((4 - math.pi) / 2)

VELOCITY_ORDER = 1
ACCELERATION_ORDER = 2
SEASONALITY_ORDER = 3


@dataclass(frozen=True)
class LeastSquaresFit:
    """Result of fitting one model to many time series.

    Attributes
    ----------
    coefficients : np.ndarray
        Shape (pixels, parameters). NaN where the fit was not possible.
    inverse_normal : np.ndarray
        ``inv(G.T @ G)`` for each pixel, shape (pixels, parameters, parameters).
    residuals : np.ndarray
        Data minus model, shape (pixels, dates). NaN for missing data.
    valid_count : np.ndarray
        Number of dates used for each pixel.
    """

    coefficients: np.ndarray
    inverse_normal: np.ndarray
    residuals: np.ndarray
    valid_count: np.ndarray

    def residual_std(self) -> np.ndarray:
        """Standard deviation of the residuals, as MATLAB's ``std`` (N - 1).

        NaN for pixels that could not be fitted.
        """
        count = np.isfinite(self.residuals).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            mean = np.nansum(self.residuals, axis=1) / count
            squared_deviation = np.nansum((self.residuals - mean[:, None]) ** 2, axis=1)
            return np.sqrt(squared_deviation / (count - 1))

    def residual_rms(self) -> np.ndarray:
        """Root mean square of the residuals. NaN for pixels not fitted."""
        count = np.isfinite(self.residuals).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.sqrt(np.nansum(self.residuals**2, axis=1) / count)


def years_between(dates: list[datetime.date], start: datetime.date) -> np.ndarray:
    """Return the time from ``start`` to each date in years.

    Parameters
    ----------
    dates : list[datetime.date]
        The dates.
    start : datetime.date
        Time zero.

    Returns
    -------
    np.ndarray
        One value per date.
    """
    return np.array([(date - start).days for date in dates]) / DAYS_PER_YEAR


def seasonal_terms_supported(dates: list[datetime.date]) -> bool:
    """Check whether the dates span enough time to estimate a seasonal term.

    Parameters
    ----------
    dates : list[datetime.date]
        All acquisition dates, sorted.

    Returns
    -------
    bool
        True if they span at least ``MINIMUM_SEASONAL_SPAN_YEARS``.
    """
    return years_between([dates[-1]], dates[0])[0] >= MINIMUM_SEASONAL_SPAN_YEARS


def build_design_matrix(
    dates: list[datetime.date], polynomial_order: int, seasonal: bool
) -> np.ndarray:
    """Build the design matrix of a polynomial plus optional annual terms.

    Parameters
    ----------
    dates : list[datetime.date]
        All acquisition dates, sorted. Time zero is the first date.
    polynomial_order : int
        Highest polynomial power.
    seasonal : bool
        Whether to add ``cos`` and ``sin`` columns with a period of one year.

    Returns
    -------
    np.ndarray
        Shape (dates, parameters). Columns: ``t**k / k!`` for k = 0 to
        ``polynomial_order``, then ``cos(2 pi tau)`` and ``sin(2 pi tau)``, where
        ``tau`` is the time since ``SEASONAL_EPOCH`` in years.
    """
    time = years_between(dates, dates[0])
    columns = [time**power / math.factorial(power) for power in range(polynomial_order + 1)]
    if seasonal:
        seasonal_time = years_between(dates, SEASONAL_EPOCH)
        columns += [np.cos(2 * np.pi * seasonal_time), np.sin(2 * np.pi * seasonal_time)]
    return np.column_stack(columns)


def fit_least_squares(design_matrix: np.ndarray, data: np.ndarray) -> LeastSquaresFit:
    """Fit one model to many time series, ignoring missing (NaN) values.

    A pixel needs at least one more valid date than the model has parameters,
    so that the residuals give a standard deviation. Otherwise its results are NaN.

    Parameters
    ----------
    design_matrix : np.ndarray
        Shape (dates, parameters).
    data : np.ndarray
        Shape (pixels, dates). NaN marks missing values.

    Returns
    -------
    LeastSquaresFit
        Coefficients, their unscaled covariance and the residuals.
    """
    n_pixels = data.shape[0]
    n_parameters = design_matrix.shape[1]
    valid = np.isfinite(data)
    valid_count = valid.sum(axis=1)

    coefficients = np.full((n_pixels, n_parameters), np.nan)
    inverse_normal = np.full((n_pixels, n_parameters, n_parameters), np.nan)

    complete = valid.all(axis=1)
    if complete.any():
        # All of these pixels share one design matrix: one matrix product solves all.
        shared_inverse = np.linalg.pinv(design_matrix.T @ design_matrix)
        coefficients[complete] = data[complete] @ (design_matrix @ shared_inverse)
        inverse_normal[complete] = shared_inverse

    with_gaps = ~complete & (valid_count > n_parameters)
    if with_gaps.any():
        weights = valid[with_gaps].astype(float)
        filled = np.where(valid[with_gaps], data[with_gaps], 0.0)
        normal = np.einsum("pt,ti,tj->pij", weights, design_matrix, design_matrix)
        right_side = np.einsum("pt,ti->pi", filled, design_matrix)
        pixel_inverse = np.linalg.pinv(normal, hermitian=True)
        coefficients[with_gaps] = np.einsum("pij,pj->pi", pixel_inverse, right_side)
        inverse_normal[with_gaps] = pixel_inverse

    too_few_dates = valid_count <= n_parameters
    coefficients[too_few_dates] = np.nan
    inverse_normal[too_few_dates] = np.nan

    residuals = np.where(valid, data - coefficients @ design_matrix.T, np.nan)
    return LeastSquaresFit(coefficients, inverse_normal, residuals, valid_count)


def _seasonal_amplitude_and_peak(
    fit: LeastSquaresFit, component_std: np.ndarray
) -> dict[str, np.ndarray]:
    """Turn the cosine and sine coefficients into amplitude and peak day.

    Parameters
    ----------
    fit : LeastSquaresFit
        A fit whose last two parameters are the cosine and sine terms.
    component_std : np.ndarray
        Standard deviation scale (EGMS uses the RMSE) for each pixel.

    Returns
    -------
    dict[str, np.ndarray]
        ``seasonality``, ``seasonality_std``, ``season_peak_doy`` and
        ``season_peak_doy_std`` for each pixel.
    """
    cosine, sine = fit.coefficients[:, -2], fit.coefficients[:, -1]
    amplitude = np.hypot(cosine, sine)
    # Average variance of the two components, as in EGMS section 11.4.
    component_variance = (fit.inverse_normal[:, -2, -2] + fit.inverse_normal[:, -1, -1]) / 2
    gaussian_std = np.sqrt(component_variance) * component_std

    # c cos(x) + s sin(x) = A cos(x - phase), which peaks where x = phase.
    phase = np.arctan2(sine, cosine)
    with np.errstate(divide="ignore", invalid="ignore"):
        phase_std = gaussian_std / amplitude
    return {
        "seasonality": amplitude,
        "seasonality_std": RAYLEIGH_STD_FACTOR * gaussian_std,
        "season_peak_doy": np.mod(phase / (2 * np.pi), 1.0) * DAYS_PER_YEAR,
        "season_peak_doy_std": phase_std / (2 * np.pi) * DAYS_PER_YEAR,
    }


def model_temporal_coherence(
    residuals_mm: np.ndarray, wavelength_m: float
) -> np.ndarray:
    """Compute the temporal coherence of model residuals, as in EGMS.

    Parameters
    ----------
    residuals_mm : np.ndarray
        Residual displacement in mm, shape (pixels, dates). NaN is ignored.
    wavelength_m : float
        Radar wavelength in meters, to convert displacement to phase.

    Returns
    -------
    np.ndarray
        ``|mean(exp(j * residual phase))|`` for each pixel, between 0 and 1.
        NaN for pixels without residuals.
    """
    meters_per_mm = 1e-3
    residual_phase = 4 * np.pi / wavelength_m * residuals_mm * meters_per_mm
    phasors = np.exp(1j * residual_phase)
    count = np.isfinite(residuals_mm).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.abs(np.nansum(phasors, axis=1)) / count


def estimate_model_parameters(
    dates: list[datetime.date],
    displacement_mm: np.ndarray,
    wavelength_m: float,
    seasonal: bool,
) -> dict[str, np.ndarray]:
    """Estimate the EGMS model parameters for many time series.

    Parameters
    ----------
    dates : list[datetime.date]
        All acquisition dates, sorted.
    displacement_mm : np.ndarray
        Displacement in mm, shape (pixels, dates). NaN marks missing values.
    wavelength_m : float
        Radar wavelength in meters, for ``model_temporal_coherence``.
    seasonal : bool
        Whether the models include annual terms. Without them, the seasonal
        parameters are NaN. See ``seasonal_terms_supported``.

    Returns
    -------
    dict[str, np.ndarray]
        One value per pixel for each model column in ``attributes.py``.
        NaN where a pixel has too few valid dates for a model.
    """
    velocity_fit = fit_least_squares(
        build_design_matrix(dates, VELOCITY_ORDER, seasonal), displacement_mm
    )
    acceleration_fit = fit_least_squares(
        build_design_matrix(dates, ACCELERATION_ORDER, seasonal), displacement_mm
    )
    seasonality_fit = fit_least_squares(
        build_design_matrix(dates, SEASONALITY_ORDER, seasonal), displacement_mm
    )

    velocity_column, acceleration_column = 1, 2
    parameters = {
        "mean_velocity": velocity_fit.coefficients[:, velocity_column],
        "mean_velocity_std": np.sqrt(
            velocity_fit.inverse_normal[:, velocity_column, velocity_column]
        ) * velocity_fit.residual_std(),
        "acceleration": acceleration_fit.coefficients[:, acceleration_column],
        "acceleration_std": np.sqrt(
            acceleration_fit.inverse_normal[:, acceleration_column, acceleration_column]
        ) * acceleration_fit.residual_std(),
        "rmse": seasonality_fit.residual_rms(),
        "model_temporal_coherence": model_temporal_coherence(
            velocity_fit.residuals, wavelength_m
        ),
    }

    if seasonal:
        parameters |= _seasonal_amplitude_and_peak(
            seasonality_fit, seasonality_fit.residual_rms()
        )
    else:
        no_value = np.full(displacement_mm.shape[0], np.nan)
        for name in ("seasonality", "seasonality_std", "season_peak_doy",
                     "season_peak_doy_std"):
            parameters[name] = no_value
    return parameters
