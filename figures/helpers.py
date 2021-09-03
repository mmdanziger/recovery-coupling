"""Functions to help plot figures for recovery coupling."""
import os.path
from typing import Optional, Tuple
import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy.signal import savgol_filter

DATA_FILENAME = "outages_util_2H_Aug12_2020.csv"


def load_dfstart(fname: Optional[str] = None) -> pd.DataFrame:
    if fname is None:
        fname = os.path.join(os.path.dirname(__file__), "..", "data", DATA_FILENAME)
    dfstart = pd.read_csv(
        fname, parse_dates=["timestamp"], index_col=["util", "timestamp"]
    )
    return dfstart

def non_zero_bounds(x, y):
    return [min(x[x > 0]), max(x), min(y[y > 0]), max(y)]

def get_weighted_errorbar_line(x, y, nbins=30):
    bins = np.logspace(0, np.log10(max(x)), nbins)
    n, _ = np.histogram(x, bins=bins)
    sy, _ = np.histogram(x, bins=bins, weights=y)
    sy2, _ = np.histogram(x, bins=bins, weights=y * y)
    mean = sy / n
    std = np.sqrt(sy2 / n - mean * mean)
    bin_x = (_[1:] + _[:-1]) / 2
    return bin_x, mean, std

def get_events_ends_from_dfstart(
    dfstart: pd.DataFrame, util: str, rolling: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    subdf = dfstart.loc[util] if util else dfstart
    if rolling:
        X = subdf.rollingevents
        Y = subdf.rollingends
    else:
        X = subdf.events
        Y = subdf.ends
    X, Y = np.array(list(zip(*sorted(zip(X.values, Y.values)))))
    nandex = np.isnan(X) | np.isnan(Y)
    X = X[~nandex]
    Y = Y[~nandex]
    return X, Y


def find_linear_until_point(
    dfstart: pd.DataFrame,
    util: str,
    outputall: bool = False,
    smoothed_max: bool = None,
    topfivemax: bool = True,
) -> Tuple[float, float]:
    if smoothed_max is None:
        smoothed_max = {"window": 21, "order": 3}
    dfxy = dfstart.loc[util][["rollingevents", "rollingends"]]
    dfxy.columns = ["X", "Y"]
    nanindex = dfxy.X.isna() | dfxy.Y.isna()

    out = []
    maxres = (0, 0)
    qvec = np.linspace(0.4, 1, 200)
    for quantile in qvec:
        Xf = dfxy[~nanindex & (dfxy.X < dfxy.X.quantile(quantile))]["X"].values
        Yf = dfxy[~nanindex & (dfxy.X < dfxy.X.quantile(quantile))]["Y"].values
        res = sm.OLS(Yf, sm.add_constant(Xf)).fit()

        out.append(res)
        maxres = max(maxres, (res.rsquared, quantile))
    y = [res.rsquared for res in out]
    if smoothed_max:
        y = savgol_filter(y, smoothed_max["window"], smoothed_max["order"])
        maxidx = np.argmax(y)
        maxres = (out[maxidx].rsquared, qvec[maxidx])
    if (
        topfivemax
    ):  # if one of the almost top ones has a higher q value, take it, don't get caught on a local minima that edges above
        topfive = sorted(zip(y, qvec))[-5:]
        maxres = max(topfive, key=lambda x: x[1])
    if outputall:
        return qvec, out, maxres
    return maxres


def get_Y_Yhat(subdf: pd.DataFrame, util: str, rolling: bool, G=1):
    X, Y = get_events_ends_from_dfstart(subdf, util, rolling)
    r2, q = find_linear_until_point(subdf, util, False)
    Xf = X[X < np.quantile(X, q)].flatten()
    Yf = Y[X < np.quantile(X, q)].flatten()
    res = sm.OLS(Yf, sm.add_constant(Xf)).fit()
    Yhat = sm.add_constant(X) @ res.params.T
    Yshift = G * (Y - res.params[0]) / res.params[1]
    data_for_df = np.vstack((X, Y, Yhat, Y - Yhat, (Y) / Yhat, Yshift)).T
    df = pd.DataFrame(
        data_for_df,
        columns=["X", "Y", "Yhat", "residual", "residual_fraction", "Yshift"],
    )
    return df
