import datetime
import matplotlib
import pandas as pd
import numpy as np
from dateutil.rrule import DAILY
import matplotlib.pyplot as plt
import seaborn as sb
import matplotlib.dates as md
from pandas.plotting import register_matplotlib_converters
import helpers

register_matplotlib_converters()


sb.set("paper", font_scale=1.9)


fig_w, fig_h = (5.5, 5.5)


def plot_imelda_compared_to_weeks(util, subdf_imelda, subdf_other):
    plt.figure(figsize=(4, 4))
    xni, yni = subdf_other[["events", "ends"]].values.T
    xi, yi = subdf_imelda[["events", "ends"]].values.T
    plt.plot(xni, yni, ".", ms=6, color="limegreen")
    plt.plot(xi, yi, ".", ms=6, color="darkred")
    plt.xlabel("Outages/2H")
    plt.ylabel("Repairs/2H")
    plt.tight_layout()
    plt.savefig(f"{util}_imelda_compared_to_weeks.png")
    plt.savefig(f"{util}_imelda_compared_to_weeks.pdf")


def plot_imelda_compared_to_quantiles(util, df_imelda, df_all_not):
    df_not_clipped = df_all_not[df_all_not.events <= df_imelda.events.max()]

    df_not_clipped["event_cuts"] = pd.cut(df_not_clipped.events, 20)
    x, y = (
        df_not_clipped.groupby("event_cuts")
        .ends.quantile([0.25, 0.5, 0.75])
        .xs(0.5, level=1)
        .reset_index()
        .values.T
    )
    x, y1 = (
        df_not_clipped.groupby("event_cuts")
        .ends.quantile([0.25, 0.5, 0.75])
        .xs(0.25, level=1)
        .reset_index()
        .values.T
    )
    x, y2 = (
        df_not_clipped.groupby("event_cuts")
        .ends.quantile([0.25, 0.5, 0.75])
        .xs(0.75, level=1)
        .reset_index()
        .values.T
    )
    x = [i.mid for i in x]

    plt.figure(figsize=(4, 4))
    xi, yi = df_imelda[["events", "ends"]].values.T
    plt.fill_between(x, y1.astype(float), y2.astype(float), alpha=0.5)
    plt.plot(xi, yi, ".", ms=6, color="darkred")
    plt.xlabel("Outages/2H")
    plt.ylabel("Repairs/2H")
    plt.savefig(f"{util}_imelda_compared_to_quantiles.png")
    plt.savefig(f"{util}_imelda_compared_to_quantiles.pdf")


def plot_imelda_ends_hist(util, subdf_not, subdf_imelda):
    plt.figure(figsize=(5, 5))
    maxbin = np.log10(subdf_not.ends.max())
    bins = np.logspace(1.1, maxbin, 35)
    subdf_not.ends.hist(bins=bins, density=True, label="Not Imelda", alpha=0.5)
    subdf_imelda.ends.hist(bins=bins, density=True, label="Imelda", alpha=0.5)
    plt.legend(loc="best")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Repairs executed in 2H window")
    plt.ylabel("Normalized Frequency")
    plt.tight_layout()
    plt.savefig(f"{util}_imelda_ends_hist.png")
    plt.savefig(f"{util}_imelda_ends_hist.pdf")


def split_imelda_and_comparison(df):
    df.index.name = "timestamp"
    df_all_not = df.reset_index()[
        (df.reset_index().timestamp < "2019-09-18")
        | ((df.reset_index().timestamp > "2019-09-25"))
    ].set_index("timestamp")
    df_imelda = df["2019-09-18":"2019-09-25"]
    df_near_not = pd.concat(
        [df["2019-09-11":"2019-09-18"], df["2019-09-25":"2019-10-02"]]
    )
    return df_imelda, df_all_not, df_near_not


def plot_imelda_elastic_response(dfstart, util):
    plt.figure(figsize=(4, 4))
    qvec, out, [r, q] = helpers.find_linear_until_point(dfstart, util, True)
    dfxy = helpers.plot_left_linear(dfstart, util, quantile=q, scatter=False)
    plt.scatter(
        dfxy.loc["2019-09-18":"2019-09-22"].X.values.flatten(),
        dfxy.loc["2019-09-18":"2019-09-22"].Y.values.flatten(),
        color="darkred",
    )
    xp, yp, score = helpers.spline_fit(
        *helpers.get_events_ends_from_dfstart(dfstart, util, True), -1
    )
    plt.plot(xp, yp)
    plt.axis([0, 2000, 0, 600])
    plt.tight_layout()

    plt.savefig(f"{util}_imelda_highlighted_response.png")
    plt.savefig(f"{util}_imelda_highlighted_response.pdf")


def get_day_of_week_hourly_average_outages(df_outages, storm_time_sequence):
    """Gets average outages for the day of week and hour of day from non-imelda period.

    Args:
        df_outages (pd.DataFrame): dataframe of outages
        storm_time_sequence (list of datetime.datetime): sequence of storm times

    Returns:
        np.ndarray, np.ndarray, np.ndarray: 0.75, 0.5 and 0.25 quantile outages per hour,
            aligned with the sequence of times provided in `storm_time_sequence`
    """
    dfstart1H = helpers.get_outages_per_hour(df_outages, "1H")
    dfe = (
        dfstart1H.loc[["centerpoint_tx", "entergy_latxak"]]
        .groupby(level=1)
        .sum()
        .events
    )
    dfe.index.name = "timestamp"
    dfe = dfe.reset_index()
    dfe = dfe[(dfe.timestamp < "2019-09-17") | (dfe.timestamp > "2019-09-24")]
    whfxn = lambda x: x.weekday() * 24 + x.hour
    dfe["wh"] = dfe.timestamp.apply(whfxn)
    dfwh = dfe.groupby("wh").quantile([0.25, 0.5, 0.75]).unstack()
    dfwh.columns = ["q1", "q2", "q3"]
    dictwh = dfwh.to_dict()

    def pull_column(key):
        return lambda i: dictwh[key][whfxn(i)] if whfxn(i) in dictwh[key] else np.nan

    y1 = np.vectorize(pull_column("q3"))(storm_time_sequence)
    y0 = np.vectorize(pull_column("q1"))(storm_time_sequence)
    y1 = np.ma.masked_invalid(y1)
    y0 = np.ma.masked_invalid(y0)
    ymed = np.vectorize(pull_column("q2"))(storm_time_sequence)
    return y1, ymed, y0


def plot_imelda_timeseries(
    df_outages, dftx, figure_start_ts, figure_end_ts, road_line_color="red"
):
    df_storm_outages = df_outages[
        (df_outages.start_timestamp > figure_start_ts)
        & (df_outages.end_timestamp < figure_end_ts)
    ]
    dftxtag = dftx[
        (dftx.start_timestamp > figure_start_ts)
        & (dftx.start_timestamp < figure_end_ts)
    ]

    outage_times, outage_events = helpers.get_full_event_sequence(
        df_storm_outages.loc[["centerpoint_tx", "entergy_latxak"]]
    )
    outage_datetimes = [datetime.datetime.fromtimestamp(_) for _ in outage_times]

    tr, er = helpers.get_full_event_sequence(dftxtag)

    (
        baseline_q75,
        baseline_median,
        baseline_q25,
    ) = get_day_of_week_hourly_average_outages(df_outages, outage_datetimes)

    plt.figure(figsize=(14, 4))
    plt.plot(outage_datetimes, outage_events, color="orange")
    plt.fill_between(
        outage_datetimes, baseline_q25, baseline_q75, alpha=0.3, color="orange"
    )
    plt.plot(outage_datetimes, baseline_median, alpha=0.3)
    ax = plt.gca()
    ax.yaxis.label.set_color("orange")
    plt.xticks(rotation=30)
    ax.xaxis_date()
    ax.autoscale()
    ax.grid(False)
    plt.ylabel("Outages")
    plt.twinx()
    plt.xticks(rotation=30)
    plt.plot(
        [datetime.datetime.fromtimestamp(_) for _ in tr], er, color=road_line_color
    )
    plt.ylabel("Flooded Roads")
    ax = plt.gca()
    ax.yaxis.label.set_color(road_line_color)
    ax.xaxis_date()
    ax.autoscale()
    ax.grid(False)
    loc = md.AutoDateLocator(minticks=10)
    loc.intervald[DAILY] = [1]
    fmt = md.AutoDateFormatter(loc)
    fmt = md.DateFormatter("%a %m-%d")
    ax.xaxis.set_major_locator(loc)
    ax.xaxis.set_major_formatter(fmt)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30)
    xmin = datetime.datetime(2019, 9, 18).toordinal()
    xmax = datetime.datetime(2019, 9, 25).toordinal()
    if matplotlib.__version__ > "3.3.0":
        new_ordinal_correction = md.date2num(np.datetime64("0000-12-31"))
        xmin += new_ordinal_correction
        xmax += new_ordinal_correction
    plt.xlim(xmin=xmin, xmax=+xmax)
    plt.tight_layout()
    plt.savefig(f"imelda_timeseries.png")
    plt.savefig(f"imelda_timeseries.pdf")


if __name__ == "__main__":
    imelda_outages = helpers.load_imelda_power_df()
    imelda_roads = helpers.load_tx_511_imelda_df(aggregate_to_patch=False)

    plot_imelda_timeseries(
        imelda_outages,
        imelda_roads,
        figure_start_ts=datetime.datetime(2019, 9, 10).timestamp(),
        figure_end_ts=datetime.datetime(2019, 10, 1).timestamp(),
        road_line_color="red",
    )
    imelda_outages_agg = helpers.get_outages_per_hour(imelda_outages, "2H")

    for util in ["entergy_latxak", "centerpoint_tx"]:
        imelda_util, nonimelda_util, nonimelda_compare = split_imelda_and_comparison(
            imelda_outages_agg.loc[util]
        )

        plot_imelda_elastic_response(imelda_outages_agg, util)
        plot_imelda_compared_to_weeks(util, imelda_util, nonimelda_compare)
        plot_imelda_ends_hist(util, nonimelda_util, imelda_util)
