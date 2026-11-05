"""Tests for the EGMS-style model fits, on synthetic time series."""

import datetime
import math

import numpy as np
import pytest

from postprocessing.timeseries_model import (
    DAYS_PER_YEAR,
    SEASONAL_EPOCH,
    build_design_matrix,
    estimate_model_parameters,
    fit_least_squares,
    seasonal_terms_supported,
)

WAVELENGTH_M = 0.0555
# Three years of 12-day acquisitions.
DATES = [datetime.date(2021, 1, 3) + datetime.timedelta(days=12 * i) for i in range(92)]
VELOCITY = 5.0  # mm/year
ACCELERATION = 2.0  # mm/year^2
AMPLITUDE = 4.0  # mm
PEAK_DAY = 200.0  # day of the year


def _synthetic_series(dates, noise_mm=0.0, seed=0):
    """Create one series with known velocity, acceleration and seasonal term."""
    time = np.array([(d - dates[0]).days for d in dates]) / DAYS_PER_YEAR
    calendar_time = np.array([(d - SEASONAL_EPOCH).days for d in dates]) / DAYS_PER_YEAR
    series = (
        VELOCITY * time
        + 0.5 * ACCELERATION * time**2
        + AMPLITUDE * np.cos(2 * np.pi * (calendar_time - PEAK_DAY / DAYS_PER_YEAR))
    )
    noise = np.random.default_rng(seed).normal(0, noise_mm, len(dates))
    return series + noise


def test_noise_free_series_returns_the_true_parameters():
    parameters = estimate_model_parameters(
        DATES, _synthetic_series(DATES)[None, :], WAVELENGTH_M, seasonal=True
    )
    assert parameters["acceleration"][0] == pytest.approx(ACCELERATION, abs=1e-6)
    assert parameters["seasonality"][0] == pytest.approx(AMPLITUDE, abs=1e-6)
    assert parameters["season_peak_doy"][0] == pytest.approx(PEAK_DAY, abs=1e-3)
    assert parameters["rmse"][0] == pytest.approx(0, abs=1e-6)
    assert parameters["model_temporal_coherence"][0] < 1  # acceleration not in model


def test_mean_velocity_is_the_linear_term_of_linear_plus_seasonal_fit():
    """EGMS's mean velocity absorbs the acceleration as the best straight line."""
    series = _synthetic_series(DATES)
    design = build_design_matrix(DATES, polynomial_order=1, seasonal=True)
    expected = np.linalg.lstsq(design, series, rcond=None)[0][1]
    parameters = estimate_model_parameters(DATES, series[None, :], WAVELENGTH_M, True)
    assert parameters["mean_velocity"][0] == pytest.approx(expected)


def test_standard_deviations_follow_the_egms_formulas():
    series = _synthetic_series(DATES, noise_mm=2.0)
    parameters = estimate_model_parameters(DATES, series[None, :], WAVELENGTH_M, True)

    design = build_design_matrix(DATES, polynomial_order=1, seasonal=True)
    inverse = np.linalg.inv(design.T @ design)
    residuals = series - design @ (inverse @ design.T @ series)
    expected_velocity_std = math.sqrt(inverse[1, 1]) * np.std(residuals, ddof=1)
    assert parameters["mean_velocity_std"][0] == pytest.approx(expected_velocity_std)

    design = build_design_matrix(DATES, polynomial_order=3, seasonal=True)
    inverse = np.linalg.inv(design.T @ design)
    residuals = series - design @ (inverse @ design.T @ series)
    rmse = math.sqrt(np.mean(residuals**2))
    expected_amplitude_std = (
        math.sqrt((4 - math.pi) / 2 * (inverse[4, 4] + inverse[5, 5]) / 2) * rmse
    )
    assert parameters["rmse"][0] == pytest.approx(rmse)
    assert parameters["seasonality_std"][0] == pytest.approx(expected_amplitude_std)


def test_gaps_give_the_same_result_as_dropping_the_dates():
    series = _synthetic_series(DATES, noise_mm=1.0)
    with_gaps = series.copy()
    with_gaps[[5, 40, 41]] = np.nan
    kept = [i for i in range(len(DATES)) if i not in (5, 40, 41)]

    design = build_design_matrix(DATES, polynomial_order=2, seasonal=True)
    both = fit_least_squares(design, np.vstack([series, with_gaps]))
    dropped = fit_least_squares(design[kept], series[kept][None, :])
    np.testing.assert_allclose(both.coefficients[1], dropped.coefficients[0])
    assert both.valid_count.tolist() == [len(DATES), len(kept)]


def test_too_few_dates_give_nan():
    series = np.full((1, len(DATES)), np.nan)
    series[0, :4] = [0.0, 1.0, 2.0, 3.0]  # linear + seasonal has 4 parameters
    parameters = estimate_model_parameters(DATES, series, WAVELENGTH_M, seasonal=True)
    assert np.isnan(parameters["mean_velocity"][0])
    assert np.isnan(parameters["mean_velocity_std"][0])
    assert np.isnan(parameters["model_temporal_coherence"][0])


def test_without_seasonal_terms_seasonal_columns_are_nan():
    short = DATES[:10]
    parameters = estimate_model_parameters(
        short, _synthetic_series(short)[None, :], WAVELENGTH_M, seasonal=False
    )
    assert np.isnan(parameters["seasonality"][0])
    assert np.isfinite(parameters["mean_velocity"][0])


def test_seasonal_terms_need_one_year():
    assert seasonal_terms_supported(DATES)
    assert not seasonal_terms_supported(DATES[:20])


def test_perfect_linear_series_has_model_coherence_one():
    linear = np.array([(d - DATES[0]).days for d in DATES]) / DAYS_PER_YEAR * VELOCITY
    parameters = estimate_model_parameters(DATES, linear[None, :], WAVELENGTH_M, True)
    assert parameters["model_temporal_coherence"][0] == pytest.approx(1.0)
    assert parameters["mean_velocity"][0] == pytest.approx(VELOCITY)
