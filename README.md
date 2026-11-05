# MintPy post-processing: one point file per run

Turns the output of an [isce2-dolphin-mintpy](../isce2-dolphin-minpty) run
(ISCE2 coregistration → dolphin phase linking → MintPy corrections) into **one
GeoParquet file**, with one point per radar pixel. Each point holds:

- its exact position (from the ISCE2 geometry), viewing geometry and water fraction,
- MintPy's and dolphin's quality measures,
- the displacement time series (`D<YYYYMMDD>` columns, in mm), after all MintPy
  corrections (troposphere, DEM error, ramps, ... as configured),
- time-series model parameters computed **as defined by EGMS**: mean velocity,
  acceleration and seasonal amplitude with standard deviations, RMSE and the day
  of the year of the seasonal peak.

It is the counterpart of [dolphin_postprocessing](../dolphin_postprocessing) for
the compass-dolphin pipeline. Column names, the file name pattern, the metadata
layout and the time-series model are the same where possible, so results of both
pipelines can be compared directly.

## Setup

```bash
mamba env create -f environment.yaml
mamba activate mintpy-postprocessing
```

This environment does not contain ISCE2, dolphin or MintPy. The export only reads
their HDF5 and GeoTIFF output.

## How to run

```bash
python run_postprocessing.py --run-dir ~/projects/ADUCAT/my_repos/isce2-dolphin-minpty/runs/test_asc073 --area-name vienna
```

`--run-dir` is the pipeline's `paths.work_dir`, the folder that contains `mintpy/`,
`dolphin/`, `isce/`, `logs/` and `pipeline_config.yaml`. The export writes two files
there, named after their content:

```
S1_ASC_073_IW1_VV_20220101_20220125_vienna.parquet         the points
S1_ASC_073_IW1_VV_20220101_20220125_vienna.metadata.json   readable copy of the metadata
```

Options:

| Option | Meaning |
| --- | --- |
| `--run-dir` | The pipeline's work directory |
| `--timeseries-file` | MintPy time series to export. Default: the end of MintPy's correction chain (see below) |
| `--output-dir` | Folder for the output, which is named automatically. Default: `--run-dir` |
| `--output` | Output file, if you want to choose the name yourself |
| `--area-name` | Name of the area, e.g. `vienna`. Added to the file name and metadata |
| `--exclude NAME ...` | Groups or columns to leave out (see below) |
| `--all-pixels` | Write every pixel, not only those in MintPy's final mask |
| `--wavelength` | Radar wavelength in m. Default: from the time series |
| `--block-rows` | Radar rows read at a time (default 256). Smaller uses less memory |
| `--overwrite` | Replace an existing output file |

Read the result with `geopandas.read_parquet(path)`, or open it in QGIS.

### Which pixels are written

By default, only pixels in **MintPy's final mask** (`maskTempCoh.h5`) are written.
That is MintPy's counterpart of dolphin's own masking. In this pipeline, the mask
contains pixels that pass both MintPy's temporal coherence threshold and dolphin's
phase-linking coherence threshold, and that lie in a valid connected component.
Pixels masked as water before unwrapping fail the connected component test.

`--all-pixels` writes every pixel with a valid position. The mask is still
written as the column `mask_temporal_coherence`, so you can filter later.

### Which time series is written

MintPy writes one time series file per correction step, and each name extends the
previous one: `timeseries.h5` → `timeseries_ERA5.h5` → `timeseries_ERA5_demErr.h5`.
The export uses the end of that chain, which is the fully corrected series that MintPy's
`velocity.h5` is also computed from. If the files don't form one chain, for example
after rerunning MintPy with another tropospheric correction in the same folder,
choose one with `--timeseries-file`.

## Run with Docker

```bash
docker build -t mintpy-postprocessing .

docker run --rm --user "$(id -u):$(id -g)" \
  -v ~/projects/ADUCAT/my_repos/isce2-dolphin-minpty/runs/test_asc073:/data/run \
  ghcr.io/lukas-scharf-plus/mintpy_postprocessing:latest \
  --run-dir /data/run --area-name vienna
```

The image's entry point is `run_postprocessing.py`, so options go straight after the
image name. On every push to `main`, GitHub Actions (`.github/workflows/docker.yml`)
builds the image, runs the tests in it and publishes it to
`ghcr.io/lukas-scharf-plus/mintpy_postprocessing` with the tags `latest` and
`sha-<commit>`. A git tag `v1.0.0` also publishes `1.0.0` and `1.0`. The commit is
stored in every output file's metadata.

Run the tests inside the image:

```bash
docker run --rm --entrypoint /usr/local/bin/_entrypoint.sh \
  mintpy-postprocessing python -m pytest /app/tests
```

## How to run the tests

```bash
python -m pytest tests
```

The tests build a tiny synthetic pipeline run and need no real data.

## Output columns

Every column stores its unit and description in the Parquet file. `point_id` and
the point geometry are always written. Points are in WGS84 longitude/latitude
(`OGC:CRS84`).

| Group (`--exclude`) | Columns | Source |
| --- | --- | --- |
| `position` | `point_id`, `row`, `col`, `lon`, `lat` | radar grid, `mintpy/inputs/geometryRadar.h5` |
| `geometry` | `height`, `incidence_angle`, `azimuth_angle`, `los_east`, `los_north`, `los_up`, `shadow_layover` | `geometryRadar.h5` (`los_*` computed) |
| `water` | `water_fraction`, `is_water` | `geometryRadar.h5` (only with a water mask in the pipeline) |
| `mintpy_quality` | `temporal_coherence`, `avg_spatial_coherence`, `mask_temporal_coherence`, `conncomp_ok`, `num_inverted_ifgrams` | `temporalCoherence.h5`, `avgSpatialCoh.h5`, `maskTempCoh.h5`, `maskConnComp.h5`, `numInvIfgram.h5` |
| `phase_linking` | `phase_linking_coherence`, `similarity`, `shp_count` | `dolphinTemporalCoherence.h5`, `dolphin/interferograms/` |
| `persistent_scatterers` | `is_ps`, `amp_dispersion` | `dolphin/interferograms/ps_mask_looked.tif`, `amp_dispersion_looked.tif` |
| `mintpy_velocity` | `velocity_mintpy`, `velocity_mintpy_std` (mm/year) | `velocity.h5` |
| `dem_error` | `dem_error` (m) | `demErr.h5` (only with MintPy's topographic residual correction) |
| `model` | `mean_velocity`, `mean_velocity_std`, `acceleration`, `acceleration_std`, `seasonality`, `seasonality_std`, `season_peak_doy`, `season_peak_doy_std`, `rmse` | fitted, see below |
| `model_coherence` | `model_temporal_coherence` | fitted, see below |
| `displacement` | `D<YYYYMMDD>` (mm, positive toward the satellite) | MintPy's final time series |
| `displacement_std` | `D<YYYYMMDD>_std` (mm) | `dolphin/interferograms/crlb_<date>.tif` |

A column whose source file is missing is left out with a warning. You can also exclude
single columns by name, for example `--exclude similarity`.

Notes on some columns:

- **Radar grid.** MintPy and dolphin work on the same radar grid, pixel for pixel.
  Every pixel gets its own latitude and longitude from ISCE2, so nothing is resampled.
  MintPy's geocoded `geo/` files are not used. `row` and `col` refer to the
  AOI-cropped, strided radar grid of the run.
- **`D<date>`** are relative to MintPy's reference date (`D<ref> = 0`) **and to its
  reference point** (`product.reference_point` in the metadata).
- **`temporal_coherence`** is MintPy's network-inversion coherence. Phase-linked
  interferograms have almost no closure phase, so it is close to 1 nearly everywhere.
  `phase_linking_coherence`, dolphin's temporal coherence, is the more useful quality
  measure. When the pipeline ran with `mintpy.inversion: dolphin`, both are dolphin's.
- **`velocity_mintpy`** is MintPy's linear fit of the final time series. `mean_velocity`
  adds seasonal terms when the series spans more than a year, so the two differ slightly.
- **`D<date>_std`** comes from dolphin's Cramér-Rao lower bound (CRLB) and covers phase
  linking only. Unwrapping, atmospheric and correction errors are not included.
- **`height`** is above the WGS84 ellipsoid, because the pipeline converts the DEM.
  It is the DEM height, not corrected by `dem_error`.
- **`los_east`, `los_north`, `los_up`** form the unit vector from the ground to the
  satellite. They follow MintPy's convention, where `azimuth_angle` is measured from
  north, anti-clockwise. A displacement `d` in east/north/up projects to the
  line of sight as `d_e·los_east + d_n·los_north + d_u·los_up`.

## Output file name and metadata

The file name follows the dolphin_postprocessing pattern:

```
S1_ASC_073_IW1_VV_20220101_20220125_vienna.parquet
   │   │   │   │  │        │        └ --area-name (optional)
   │   │   │   │  │        └ last acquisition
   │   │   │   │  └ first acquisition
   │   │   │   └ polarization
   │   │   └ processed subswaths
   │   └ track (relative orbit)
   └ orbit direction
```

The metadata is stored as JSON under the key `insar_processing` in the Parquet file,
and next to it as `<name>.metadata.json`:

```python
import json
import pyarrow.parquet as pq

metadata = json.loads(pq.read_schema(path).metadata[b"insar_processing"])
```

| Section | Content | Source |
| --- | --- | --- |
| `product` | sensor, missions, track, direction, swaths, polarization, area, all dates, reference date, **reference point** (pixel and lat/lon) | MintPy time series attributes, ISCE2 logs |
| `software` | ISCE2, dolphin, MintPy and pipeline versions; this code's commit and libraries | `logs/versions.json` (older runs: dolphin and MintPy logs) |
| `inputs` | every SLC with its satellite, orbit file and orbit type (`POEORB`/`RESORB`), count per orbit type, DEM file and its vertical datum, water mask file | ISCE2 unpack logs, `pipeline_config.yaml` |
| `settings` | the complete pipeline config, the complete dolphin config, the MintPy template, the time series file used, the radar grid, the export options | run directory, command line |
| `run` | processing time (UTC) and every environment variable starting with `INSAR_RUN_` | system clock, environment |

`inputs.orbit_types` shows at a glance whether any date used a restituted orbit
(`RESORB`). The SLC and orbit files themselves are not needed. ISCE2 logs which orbit
file it used for each SAFE, and that is read instead.

## How it works

1. **Find the files** (`postprocessing/run_outputs.py`). Only the MintPy geometry
   and time series are required. Every other file is optional.
2. **Collect the metadata** (`postprocessing/metadata.py`) and build the file name.
3. **Block loop** (`postprocessing/export.py`). The radar grid is read in blocks of
   rows. The selected pixels become points, and the time-series model is fitted in
   the same loop.
4. **Writing.** The file is written to `<output>.partial` and renamed only after
   success, so a crash never leaves a file that looks complete.

### Time-series model (`postprocessing/timeseries_model.py`)

The code is identical to dolphin_postprocessing, so the parameters of both pipelines can be compared directly.

The definitions follow the EGMS *Product Description and Format Specification*
(sections 5.1 and 11.4).

#### The idea in short

A displacement time series is a list of numbers, one per date. To compare points,
or to compare with EGMS, we describe each series with a few numbers:

- How fast does the point move on average? This is the **velocity**.
- Is the motion speeding up or slowing down? This is the **acceleration**.
- Does the point move back and forth every year, for example because buildings
  expand in summer? This is the **seasonal amplitude**, plus the time of year of
  the peak.
- How much of the series do these numbers *not* explain? This is the **RMSE**.

Each number comes from fitting a smooth curve (a "model") to the series with least
squares, which means choosing the curve that minimizes the sum of squared
differences to the data. Every number also gets a **standard deviation**, which
says how precisely the data determine it.

#### Time variables

For acquisition $i$ (date $D_i$, first date $D_0$):

$$
t_i = \frac{D_i - D_0}{365.25\ \text{days}} \qquad\text{(years since the first date, used for the polynomial)}
$$

$$
\tau_i = \frac{D_i - \text{2000-01-01}}{365.25\ \text{days}} \qquad\text{(calendar time, used for the seasonal term)}
$$

Both are "time in years", only measured from different starting points:

- $t$ starts at **the first acquisition of the run**. With this choice, $a_0$ is
  the displacement at the first date. It also keeps the numbers small, which
  keeps the $t^2$ and $t^3$ terms numerically well behaved.
- $\tau$ starts at **a fixed calendar date** (1 January 2000). Only the fraction
  after the decimal point matters for $\cos(2\pi\tau)$ and $\sin(2\pi\tau)$, and
  that fraction is how far into the calendar year the date lies (0.0 = 1 January,
  0.5 ≈ 1 July). So the seasonal peak is a real day of the year that you can
  compare between runs. If $t$ were used here, "peak at 0.3" would mean 0.3
  years after whichever date the run happened to start.

Example: a run starts on 2021-03-15. For the acquisition on 2022-07-01,
$t = 473 / 365.25 = 1.295$ (1.3 years into the run) and
$\tau = 8217 / 365.25 = 22.497$ (22 years after 2000, about halfway through the year).

The seasonal term is an annual oscillation:

$$
S(\tau) = c\,\cos(2\pi\tau) + s\,\sin(2\pi\tau) = A\,\cos\!\big(2\pi(\tau - \tau_\text{peak})\big),
\qquad A = \sqrt{c^2 + s^2},\quad \tau_\text{peak} = \frac{\mathrm{atan2}(s, c)}{2\pi}
$$

The sine/cosine form keeps the model linear in its unknowns, so least squares
solves it directly. Amplitude and peak time follow from $c$ and $s$ afterwards.

#### The three models

EGMS fits **three separate models** to each series $d(t)$ (in mm). Each parameter
comes from the simplest model that contains it:

| Model | Formula | Gives |
| --- | --- | --- |
| 1: linear + seasonal | $d(t) = a_0 + v\,t + S(\tau)$ | `mean_velocity` $= v$ |
| 2: quadratic + seasonal | $d(t) = a_0 + v\,t + \tfrac{1}{2}\,\text{acc}\,t^2 + S(\tau)$ | `acceleration` $= \text{acc}$ |
| 3: cubic + seasonal | $d(t) = a_0 + a_1 t + \tfrac{1}{2}a_2 t^2 + \tfrac{1}{6}a_3 t^3 + S(\tau)$ | `seasonality` $= A$, `season_peak_doy`, `rmse` |

Why three models and not one?

- **Velocity** comes from model 1, without acceleration. It is then the best
  straight line through the whole period: the average rate. In model 2, $v$ would
  be the rate at the first date only.
- **Acceleration** needs at least a quadratic term, so it comes from model 2.
- **Seasonality and RMSE** come from the most flexible model 3. A slow non-linear
  trend that a straight line would miss is absorbed by the cubic, so it neither
  leaks into the seasonal amplitude nor inflates the RMSE.

EGMS writes the quadratic term as $a_2 x^2$ and reports $2a_2$. Writing it as
$\tfrac{1}{2}\,\text{acc}\,t^2$ gives the same number directly.

#### Least squares and standard deviations

Each model is written as $\mathbf{d} = G\,\mathbf{x} + \mathbf{r}$, where $G$ has
one row per date and one column per term ($1$, $t$, $t^2/2$, …, $\cos$, $\sin$).
The estimate and its unscaled covariance are:

$$
\hat{\mathbf{x}} = (G^\mathsf{T}G)^{-1}G^\mathsf{T}\mathbf{d}, \qquad Q = (G^\mathsf{T}G)^{-1},
\qquad \mathbf{r} = \mathbf{d} - G\hat{\mathbf{x}}
$$

$Q$ only depends on the dates. It is scaled with the noise level, estimated from
the residuals $r_i$ of the $n$ dates:

$$
s_r = \sqrt{\frac{1}{n-1}\sum_i (r_i - \bar r)^2}, \qquad
\text{rmse} = \sqrt{\frac{1}{n}\sum_i r_i^2}\ \ \text{(model 3)}
$$

| Column | Formula |
| --- | --- |
| `mean_velocity_std` | $\sigma_v = \sqrt{Q_{vv}}\; s_r$ (model 1) |
| `acceleration_std` | $\sigma_\text{acc} = \sqrt{Q_\text{acc,acc}}\; s_r$ (model 2) |
| `seasonality_std` | $\sigma_A = \sqrt{\tfrac{4-\pi}{2}}\;\sqrt{\tfrac{Q_{cc} + Q_{ss}}{2}}\;\text{rmse}$ (model 3) |
| `season_peak_doy` | $(\tau_\text{peak} \bmod 1) \cdot 365.25$ days (0 = 1 January) |
| `season_peak_doy_std` | $\dfrac{1}{2\pi}\,\dfrac{\sqrt{(Q_{cc} + Q_{ss})/2}\;\text{rmse}}{A} \cdot 365.25$ days |

What the formulas mean:

- **More dates and a longer period make $Q$ smaller**, so the parameters become
  more precise. Noisier data (larger $s_r$) makes them less precise.
- **The factor $\sqrt{(4-\pi)/2} \approx 0.66$** appears because an amplitude can
  never be negative. Even pure noise gives some $A > 0$, and such an amplitude
  follows a Rayleigh distribution rather than a normal one.
- **The peak day becomes meaningless when the amplitude is small**, because a
  tiny oscillation has no clear peak. Its standard deviation grows as $1/A$.
  EGMS reports no phase; `season_peak_doy` is an addition here.
- **These are formal standard deviations.** They assume independent noise on each
  date. Atmospheric effects are correlated in time, so the true uncertainty is
  usually larger. EGMS notes the same ("without considering the atmospheric phase
  screen").

#### Model temporal coherence

`model_temporal_coherence` turns the residuals of model 1 back into phase and
measures how consistent they are:

$$
\gamma_\text{model} = \left|\frac{1}{n}\sum_i \exp\!\left(j\,\frac{4\pi}{\lambda}\,r_i\right)\right|
$$

It is 1 when model 1 explains the series perfectly. It drops when the residuals
are large compared with the wavelength $\lambda$, either because the data are noisy
or because the point moves in a way the model does not describe, such as a sudden
step. dolphin's `phase_linking_coherence` uses no motion model: it checks whether one
phase per date explains all interferograms. So use `phase_linking_coherence` to judge
data quality, and `model_temporal_coherence` to judge whether model 1 fits.

#### Reading the values

- Signs: displacement is positive **toward the satellite**. A negative
  `mean_velocity` means moving away from the satellite, for example subsidence.
- A parameter is clearly present when it is larger than about two to three times
  its standard deviation, for example `abs(acceleration) > 3 * acceleration_std`.
- `season_peak_doy` around 200–250 (July to early September) with a clear
  amplitude is typical for thermal expansion of buildings.

Limits:

- With **less than one year** of data, the seasonal terms cannot be separated from
  the polynomial. All models are then fitted without them, and the seasonal
  columns are empty (a warning says so).
- A pixel needs at least one more valid date than the model has parameters
  (5, 6 and 7 dates). Otherwise its model columns are empty.
- Pixels without gaps share one design matrix and are solved together with one
  matrix product. Pixels with gaps are solved separately, still vectorized.

## Project structure

- `run_postprocessing.py`: command-line entry point.
- `postprocessing/run_outputs.py`: finds the MintPy, dolphin and ISCE2 files of a run.
- `postprocessing/attributes.py`: column names, groups, units and `--exclude`.
- `postprocessing/timeseries_model.py`: EGMS-style least-squares fits (same as dolphin_postprocessing).
- `postprocessing/export.py`: block loop and GeoParquet writer.
- `postprocessing/metadata.py`: processing metadata and the output file name.
- `tests/`: automated tests. `conftest.py` builds a tiny synthetic pipeline run.
- `environment.yaml`: conda environment definition.
- `Dockerfile`, `.dockerignore`: container image built from `environment.yaml`.
- `.github/workflows/docker.yml`: builds, tests and publishes the image to ghcr.io.
- `CLAUDE.md`: working instructions for Claude Code.
