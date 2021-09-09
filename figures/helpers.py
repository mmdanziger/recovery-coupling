"""Functions to help plot figures for recovery coupling."""
import os.path
import datetime
from typing import Optional, Tuple
import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy.signal import savgol_filter
import mercantile
from matplotlib import colors
from haversine import haversine,Unit
import tqdm
from lifelines import KaplanMeierFitter

DATA_FILENAME = "outages_util_2H_Aug12_2020.csv"
IMELDA_ROADS_FILENAME = "tx_511.csv"
IMELDA_COUNTIES_OUTAGES = "imelda_counties_outages_to_Jul06_2020.csv"

def load_dfstart(fname: Optional[str] = None) -> pd.DataFrame:
    if fname is None:
        fname = os.path.join(os.path.dirname(__file__), "..", "data", DATA_FILENAME)
    dfstart = pd.read_csv(
        fname, parse_dates=["timestamp"], index_col=["util", "timestamp"]
    )
    return dfstart


def load_tx_511_imelda_df(fname: Optional[str] = None) -> pd.DataFrame:
    if fname is None:
        fname = os.path.join(
            os.path.dirname(__file__), "..", "data", IMELDA_ROADS_FILENAME
        )
    df_roads = pd.read_csv(fname)
    df_roads_tag = df_roads[
        df_roads["7"].apply(
            lambda x: any([i in x.lower() for i in ["water", "flood", "damage"]])
        )
    ]
    df_roads_tag = add_spatial_index_column(df_roads_tag)
    aggdict = {
        "start_timestamp": "min",
        "end_timestamp": "max",
        "lat": "mean",
        "lng": "mean",
    }
    df_roads_tag.groupby("spatial_index").agg(aggdict)
    df_roads_tag = df_roads_tag.reset_index()
    return df_roads_tag

def load_imelda_power_df(fname: Optional[str] = None) -> pd.DataFrame:
    if fname is None:
        fname = os.path.join(
            os.path.dirname(__file__), "..", "data", IMELDA_COUNTIES_OUTAGES
        )
    df = pd.read_csv(fname,
                     dtype={"spatial_index": str, "end_timestamp": int, "incid": str, "block": str, "block_group": str},
                     low_memory=False)
    df = df.rename(columns={"Unnamed: 0": "util"})
    df = df.set_index(["util", "index_copy"])
    df["start_timestamp_observed"] = [int(i[1].split("_")[0]) for i in df.index]
    df["start_timestamp_reported"] = df["start_timestamp"].values
    df = df.loc[["centerpoint_tx","entergy_latxak"]]
    df.end_timestamp = df.end_timestamp.apply(int)
    df.start_timestamp = df.start_timestamp.apply(int)
    timestring = lambda x: datetime.datetime.fromtimestamp(x).strftime("%Y-%m-%dT%H:%M:%S")
    df["start_timestring"] = df.start_timestamp.apply(timestring)
    df["end_timestring"] = df.end_timestamp.apply(timestring)
    return df

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


def get_neighboring_keys(qk, include_self=True):
    tile0 = mercantile.quadkey_to_tile(qk)
    neighbor_qk = []
    for dx in [1, 0, -1]:
        for dy in [1, 0, -1]:
            if include_self or not (dx == 0 and dy == 0):
                neighbor_qk.append(
                    mercantile.quadkey(tile0.x + dx, tile0.y + dy, tile0.z)
                )
    return neighbor_qk


def add_spatial_index_column(df, overwrite=False, zoom=19):
    if "spatial_index" not in df.columns:
        df["spatial_index"] = ""
    if not overwrite:
        df.loc[
            ((df.spatial_index == "") | (df.spatial_index.isna())) & ~df.lat.isna(),
            "spatial_index",
        ] = df.loc[
            ((df.spatial_index == "") | (df.spatial_index.isna())) & ~df.lat.isna()
        ].apply(
            lambda x: mercantile.quadkey(mercantile.tile(x["lng"], x["lat"], zoom)),
            axis=1,
        )
    else:
        df.loc[~df.lat.isna(), "spatial_index"] = df.loc[~df.lat.isna()].apply(
            lambda x: mercantile.quadkey(mercantile.tile(x["lng"], x["lat"], zoom)),
            axis=1,
        )
    return df

def add_zoom_column(df, z=11):
    has_spatial_index = ~df.spatial_index.isna() & ~(df.spatial_index == "")
    df.loc[has_spatial_index, "spatial_index_z%i" % z] = df.loc[has_spatial_index, "spatial_index"].transform(
        lambda x: x[:z] if len(x) > z else "")
    not_enough_z = (df["spatial_index_z%i" % z] == "") & ~df.lat.isna()
    df.loc[not_enough_z, "spatial_index_z%i" % z] = df.loc[not_enough_z].apply(
        lambda x: mercantile.quadkey(mercantile.tile(x["lng"], x["lat"], z)), axis=1)
    return df


def get_outages_per_hour(df, freq="H", window_size=3):
    df.loc[(df.start_timestamp > 10000000000), ["start_timestamp"]] = df[
        (df.start_timestamp > 10000000000)
    ].start_timestamp.apply(lambda x: x // 1000)

    start_df = (
        df.reset_index()
        .set_index(pd.DatetimeIndex(pd.to_datetime(df["start_timestamp"], unit="s")))
        .rename(columns={"level_0": "util", "level_1": "fidx"})
    )
    groupby_key = (
        ["util", pd.Grouper(freq=freq)]
        if "util" in start_df.columns
        else pd.Grouper(freq=freq)
    )
    start_per_hour = start_df.groupby(groupby_key).start_timestamp.count()

    end_df = (
        df.reset_index()
        .set_index(pd.DatetimeIndex(pd.to_datetime(df["end_timestamp"], unit="s")))
        .rename(columns={"level_0": "util", "level_1": "fidx"})
    )
    end_per_hour = -end_df.groupby(groupby_key).start_timestamp.count()
    med_recovery_per_hour = (
        start_df.assign(
            recovery=start_df["end_timestamp"] - start_df["start_timestamp"]
        )
        .groupby(groupby_key)
        .agg({"recovery": "median"})
    )
    cust_a_per_hour = start_df.groupby(groupby_key).cust_a.sum()
    cust_a_end_per_hour = -end_df.groupby(groupby_key).cust_a.sum()
    startend = pd.concat(
        [
            start_per_hour,
            end_per_hour,
            med_recovery_per_hour,
            cust_a_per_hour,
            cust_a_end_per_hour,
        ],
        axis=1,
        sort=True,
    )
    startend.columns = [
        "starts",
        "ends",
        "median_recovery",
        "cust_a_start",
        "cust_a_end",
    ]
    startend = startend.transform(np.nan_to_num)
    if "util" in start_df.columns:
        cumsum_df = startend.groupby(level=0).cumsum()
    else:
        cumsum_df = startend.cumsum()
    startend = startend.join(cumsum_df, rsuffix="_cumsum")
    startend["leftover_events"] = startend["starts_cumsum"] + startend["ends_cumsum"]
    startend["leftover_cust_a"] = (
        startend["cust_a_start_cumsum"] + startend["cust_a_end_cumsum"]
    )

    startend["events"] = (
        np.array([0] + [i for i in startend["leftover_events"].values][:-1])
        + startend["starts"]
    )
    startend["cust_a"] = (
        np.array([0] + [i for i in startend["leftover_cust_a"].values][:-1])
        + startend["cust_a_start"]
    )
    if startend.index.nlevels == 1:
        startend["hour"] = startend.index
    else:
        startend["hour"] = startend.index.get_level_values(1)
    startend["hour"] = startend.hour.apply(lambda x: x.time().hour)
    startend["rollingleftover_events"] = startend.leftover_events.rolling(
        window=window_size
    ).mean()
    startend["rollingevents"] = startend.events.rolling(window=window_size).mean()
    startend["rollingends"] = -startend.ends.rolling(window=window_size).mean()
    startend["rollingstarts"] = startend.starts.rolling(window=window_size).mean()
    startend["rollingrecovery"] = startend.median_recovery.rolling(
        window=window_size
    ).mean()
    x = startend.starts_cumsum.values
    xtag = [x[i] - x[i - 5] for i in range(5, len(x))]
    xtag = list(x[:5]) + xtag
    startend["recent_starts"] = xtag
    startend.ends = -startend.ends
    startend["events_z"] = (
        startend.events - startend.events.mean()
    ) / startend.events.std()
    return startend

def calc_distances(row,z,dfz):
    close = dfz[dfz["spatial_index_z%i"%z].isin(get_neighboring_keys(row["spatial_index"][:z]))]
    time_overlap = close[(close.start_timestamp < row["end_timestamp"]) & (close.end_timestamp > row["start_timestamp"])]
    out=[]
    for idx,other_row in time_overlap.iterrows():
        d=haversine(other_row[["lat","lng"]],row[["lat","lng"]],unit=Unit.KILOMETERS)
        out.append((idx,d))
    return out


def kmffit(dvec):
    kmf=KaplanMeierFitter()
    kmf.fit(dvec)
    return kmf

def add_closest_flooded_roads(df_roads, df_power_storm,     scanning_z=11):
    #Calculate the nearest flooded road for each outage
    df_power_storm = add_zoom_column(df_power_storm,scanning_z)
    distmatrix = np.zeros( (len(df_roads),len(df_power_storm)) ) + 10000
    for i in tqdm.tqdm(range(len(df_roads))):
        dists=calc_distances(df_roads.loc[i],scanning_z,df_power_storm)
        for j,d in dists:
            distmatrix[i,j] = d

    dmins = np.array([(i,distmatrix[:,i].argmin(),distmatrix[:,i].min()) for i in range(distmatrix.shape[1])])
    closest_flooded_road_dist = dmins[:,2]
    closest_flooded_road = dmins[:,1].astype(int)
    df_power_storm["closest_flooded_road_dist"] = closest_flooded_road_dist
    df_power_storm["closest_flooded_road"] = closest_flooded_road
    return df_power_storm

def add_outage_details_to_df_roads(df_roads, df_power_storm):
    dfz_gb_roadidx = df_power_storm[df_power_storm.closest_flooded_road_dist < 10000].groupby("closest_flooded_road").agg({"durationh":["median","mean","count"]})
    dfz_gb_roadidx.columns = ["outage_duration_median","outage_duration_mean","outage_count"]

    df_roads = add_spatial_index_column(df_roads)
    df_roads = df_roads.join(dfz_gb_roadidx)
    return df_roads

def add_distance_partition_by_nearest_flooded_road(df_power_storm, km_cuts):
    km_cut_labels = ["%i - %i"%(km_cuts[i],km_cuts[i+1]) for i in range(len(km_cuts) -2) ]
    km_cut_labels.append("> %i"%km_cuts[-2] )
    df_power_storm["distance_partition"] = pd.cut(df_power_storm.closest_flooded_road_dist,km_cuts, labels=km_cut_labels)
    return df_power_storm

def kmf_curves_by_distance_partition(df_power_storm):
    return {part : kmffit(durations_by_partition(df_power_storm, part)) for part in df_power_storm.distance_partition.unique()}

def durations_by_partition(df, part):
    return df[df.distance_partition==part].durationh.values

def get_divnorm():
    try:
        return colors.TwoSlopeNorm(1, vmin=0, vmax=2)
    except AttributeError:
        return colors.DivergingNorm(vmin=0, vcenter=1, vmax=2)
