#!/usr/bin/env python3
"""
Merge the four overlapping TPOSE-Vel assimilation runs into daily monthly
NetCDF files.

The four runs are 4-month assimilations staggered by two months, so
consecutive runs share a two-month overlap. In each overlap the outgoing run
is linearly ramped from weight 1 to 0 and the incoming run from 0 to 1 across
the overlapping days; outside overlaps a single run has weight 1. Where three
runs coincide (the two-day tail of a run that extends past its nominal
four-month span) the most recent run pair is used.

see Verdy, et al. (2023) in JGR Oceans for more detail on the smoothing.

Output is one file per variable per calendar month, subset to the upper 500 m
(where a depth dimension exists) and 20S-20N. See README.md.
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
from xmitgcm import open_mdsdataset

GRID_DIR   = '/data/SO6/TPOSE_diags/tpose6/grid_6/'
DIAGS_ROOT = '/data/SO3/edavenport/tpose6/diags'
OUTPUT_DIR = '/data/SO3/edavenport/tpose6/TPOSE-Vel'

# Final iteration of each 4-month run. ref_date is one day before the first
# daily mean so the time coordinate labels each record by the day it covers
# (iter 72 -> day 1).
WINDOWS = [
    dict(name='sep2012', run='run_iter14', ref_date='2012-08-31'),
    dict(name='nov2012', run='run_iter20', ref_date='2012-10-31'),
    dict(name='jan2013', run='run_iter14', ref_date='2012-12-31'),
    dict(name='mar2013', run='run_iter16', ref_date='2013-02-28'),
]

DELTA_T = 1200                          # seconds per model step (86400 s / 72)
ITERS   = list(range(72, 72 * 123, 72)) # 122 daily means per run

DEPTH_MAX        = 500.0                # m, retain cells with centre >= -DEPTH_MAX
LAT_MIN, LAT_MAX = -20.0, 20.0

# variable -> diagnostic file prefix
VARS = {
    'THETA': 'diag_state',
    'SALT':  'diag_state',
    'UVEL':  'diag_state',
    'VVEL':  'diag_state',
    'WVEL':  'diag_state',
    'ETAN':  'diag_surf',
}

COMPLEVEL = 4


def subset(ds):
    """Restrict to the upper DEPTH_MAX metres and LAT_MIN..LAT_MAX latitude."""
    for z in ('Z', 'Zl', 'Zu', 'Zp1'):
        if z in ds.dims:
            ds = ds.sel({z: ds[z][ds[z] >= -DEPTH_MAX]})
    for y in ('YC', 'YG'):
        if y in ds.dims:
            ds = ds.sel({y: ds[y][(ds[y] >= LAT_MIN) & (ds[y] <= LAT_MAX)]})
    return ds


def load_window(window, prefix):
    """Open one run/prefix lazily, subset in space, chunk one day per task."""
    ds = open_mdsdataset(
        data_dir=os.path.join(DIAGS_ROOT, window['name'], window['run']),
        grid_dir=GRID_DIR, iters=ITERS, prefix=[prefix],
        ref_date=window['ref_date'], delta_t=DELTA_T,
    )
    # Coordinates are big-endian float32 from disk; cast to native float so
    # label-based selection works.
    for c in ('XC', 'XG', 'YC', 'YG', 'Z', 'Zl', 'Zu', 'Zp1'):
        if c in ds.coords:
            ds[c] = ds[c].astype('float64')
    return subset(ds).chunk({'time': 1})


def build_weights(time_axes):
    """
    Return (all_times, W) where all_times is the sorted union of the run daily
    timestamps and W[i, j] is the weight of run i at all_times[j], with the
    columns summing to 1.
    """
    all_times = pd.DatetimeIndex(sorted(set().union(*[set(t) for t in time_axes])))
    n_w, n_t = len(time_axes), len(all_times)

    covered = np.zeros((n_w, n_t), dtype=bool)
    for i, t in enumerate(time_axes):
        covered[i, all_times.get_indexer(t)] = True

    overlaps = {(i, i + 1): time_axes[i].intersection(time_axes[i + 1])
                for i in range(n_w - 1)}

    W = np.zeros((n_w, n_t))
    for j in range(n_t):
        covering = np.where(covered[:, j])[0]
        if len(covering) == 1:
            W[covering[0], j] = 1.0
        else:
            # Most recent overlapping pair; linear ramp across their overlap.
            i_out, i_in = int(covering[-2]), int(covering[-1])
            ov   = overlaps[(i_out, i_in)]
            frac = ov.get_loc(all_times[j]) / (len(ov) - 1)
            W[i_out, j] = 1.0 - frac
            W[i_in,  j] = frac

    assert np.allclose(W.sum(axis=0), 1.0), 'weights do not sum to 1'
    return all_times, W


def blend(win_das, all_times, W, month_times):
    """Weighted sum of the run DataArrays over the given month timestamps."""
    tidx = all_times.get_indexer(month_times)
    out = None
    for i, da in enumerate(win_das):
        w = W[i, tidx]
        active = w > 0
        if not active.any():
            continue
        at, aw = month_times[active], w[active]
        da_a = da.sel(time=at)
        if np.all(aw == 1.0):
            contrib = da_a
        else:
            wda = xr.DataArray(aw.astype('float32'), coords={'time': at}, dims=['time'])
            contrib = wda * da_a
        if not active.all():                       # run absent on some days
            contrib = contrib.reindex(time=month_times, fill_value=0.0)
        out = contrib if out is None else out + contrib
    return out


def month_groups(all_times):
    """List of (label, month_times) for each calendar month in all_times."""
    groups = []
    for key in sorted({(t.year, t.month) for t in all_times}):
        mt = all_times[(all_times.year == key[0]) & (all_times.month == key[1])]
        label = pd.Timestamp(key[0], key[1], 1).strftime('%b%Y').lower()
        groups.append((label, mt))
    return groups


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    prefixes = sorted(set(VARS.values()))
    datasets = {p: [load_window(w, p) for w in WINDOWS] for p in prefixes}

    time_axes = [pd.DatetimeIndex(ds.time.values) for ds in datasets[prefixes[0]]]
    all_times, W = build_weights(time_axes)
    months = month_groups(all_times)

    print(f'Coverage: {all_times[0].date()} -> {all_times[-1].date()} '
          f'({len(all_times)} days, {len(months)} months)')

    for var, prefix in VARS.items():
        win_das  = [ds[var] for ds in datasets[prefix]]
        has_depth = any(z in win_das[0].dims for z in ('Z', 'Zl'))
        depth_tok = '0to500m_' if has_depth else ''
        for label, mtimes in months:
            fname = f'tpose_vel_{var}_{label}_daily_{depth_tok}20Sto20N.nc'
            outpath = os.path.join(OUTPUT_DIR, fname)
            if os.path.exists(outpath):
                print(f'skip {fname}')
                continue
            print(f'write {fname} ({len(mtimes)} days)', flush=True)
            da = blend(win_das, all_times, W, mtimes)
            da.attrs['blended_from'] = ', '.join(
                f"{w['name']}/{w['run']}" for w in WINDOWS)
            ds_out = da.to_dataset(name=var)
            ds_out.attrs['description'] = (
                'TPOSE-Vel assimilation runs merged with linear weighting '
                'across two-month overlaps; upper 500 m, 20S-20N.')
            enc = {var: {'zlib': True, 'complevel': COMPLEVEL, 'dtype': 'float32',
                         'chunksizes': (1,) + da.shape[1:]}}
            ds_out.to_netcdf(outpath, encoding=enc, format='NETCDF4')

    print('done')


if __name__ == '__main__':
    main()
