# TPOSE-Vel post-processing

`make_tpose_vel.py` merges the four overlapping TPOSE-Vel assimilation runs into
daily, monthly NetCDF files.

## Input

Four 4-month assimilation runs, each staggered from the previous by two months,
so consecutive runs share a two-month overlap. The final iteration of each run
is used:

| Run     | Directory                                            | Days covered            |
|---------|------------------------------------------------------|-------------------------|
| sep2012 | `/data/SO3/edavenport/tpose6/diags/sep2012/run_iter14` | 2012-09-01 – 2012-12-31 |
| nov2012 | `/data/SO3/edavenport/tpose6/diags/nov2012/run_iter20` | 2012-11-01 – 2013-03-02 |
| jan2013 | `/data/SO3/edavenport/tpose6/diags/jan2013/run_iter14` | 2013-01-01 – 2013-05-02 |
| mar2013 | `/data/SO3/edavenport/tpose6/diags/mar2013/run_iter16` | 2013-03-01 – 2013-06-30 |

Each run writes 122 daily-mean records (`diag_state` and `diag_surf`), MITgcm
iterations 72 to 8784 in steps of 72 (`delta_t = 1200` s, 72 steps per day).
Grid: `/data/SO6/TPOSE_diags/tpose6/grid_6`.

## Merging

Runs are combined by linear weighting over their overlaps:

- Outside any overlap a single run has weight 1.
- Within a two-month overlap the outgoing run is ramped linearly from weight 1
  to 0 and the incoming run from 0 to 1 across the overlapping days. The ramp
  reaches 1 and 0 at the first and last overlapping day respectively; weights
  sum to 1 on every day.
- Because the calendar months in each run are not all the same length, a run
  can extend one or two days past its nominal four-month span, so three runs
  briefly coincide (2013-03-01 to 2013-03-02). On such days the most recent run
  pair is used and the earliest run is dropped.

Merged coverage is 2012-09-01 to 2013-06-30 (10 months).

## Output

Written to `/data/SO3/edavenport/tpose6/TPOSE-Vel/`, one file per variable per
calendar month:

```
tpose_vel_<VAR>_<mon><year>_daily_0to500m_20Sto20N.nc
```

`ETAN` has no depth dimension, so its files omit the `0to500m` token.

- Variables: `THETA`, `SALT`, `UVEL`, `VVEL`, `WVEL`, `ETAN`.
- Daily means, one calendar month per file.
- Where a depth dimension exists, cells whose centre lies within the upper 500 m
  are retained (43 levels, `THETA`/`SALT`/`UVEL`/`VVEL` on `Z`, `WVEL` on `Zl`).
- Latitude restricted to 20S–20N (`YC`, or `YG` for `VVEL`).
- Full model longitude range is kept.
- Data stored as float32, zlib-compressed.

## Assumptions

- The highest-numbered `run_iter` in each run directory is the final,
  best-converged iteration.
- All runs share the tpose6 `grid_6` grid.
- Daily-mean records are labelled by the day they cover (the run `ref_date` is
  set one day before the first record so iteration 72 maps to day 1).
- A cell is "within 500 m" if its centre depth is >= -500 m.

## Running

```
conda activate tpose
python make_tpose_vel.py
```

Existing output files are skipped, so the script can be re-run to resume.

## Loading

```python
import xarray as xr
ds = xr.open_dataset(
    '/data/SO3/edavenport/tpose6/TPOSE-Vel/'
    'tpose_vel_THETA_sep2012_daily_0to500m_20Sto20N.nc')
```
