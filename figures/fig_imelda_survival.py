import datetime
import pandas as pd
import helpers
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sb
import statsmodels.api as sm

sb.set("paper", font_scale=1.9)
fig_w, fig_h = (5.5, 5.5)


def plot_simplified_survival_curve(df_power_storm, df_prev60):

    df_power_storm = helpers.add_distance_partition_by_nearest_flooded_road(
        df_power_storm, [0, 5, 30, 10000]
    )

    kmflist = helpers.kmf_curves_by_distance_partition(df_power_storm)

    f = plt.figure(figsize=(fig_w, fig_h))
    a = f.add_subplot(111)
    event_count = {}
    for z, kmf in sorted(
        kmflist.items(), key=lambda x: int(x[0].strip("<>").split()[0])
    ):
        event_count[z] = len(df_power_storm[df_power_storm.distance_partition == z])
        kmf.plot_survival_function(ax=a, label=f"{z} km ({event_count[z]} events)")

    helpers.kmffit(df_prev60.durationh.values).plot_survival_function(
        ax=a, label=f"non-Imelda ({len(df_prev60)} events)"
    )

    plt.xlim(xmin=0, xmax=36)
    plt.xlabel("Hours")
    plt.ylabel("Fraction of Outages Unrepaired")
    plt.tight_layout()
    plt.savefig("imelda_survival_simplified.png")
    plt.savefig("imelda_survival_simplified.pdf")


def plot_imelda_survival(df_power_storm):
    # partition the data
    df_power_storm = helpers.add_distance_partition_by_nearest_flooded_road(
        df_power_storm, [0, 2, 5, 10, 20, 10000]
    )
    kmflist = helpers.kmf_curves_by_distance_partition(df_power_storm)
    f = plt.figure(figsize=(fig_w, fig_h))
    a = f.add_subplot(111)
    event_count = {}
    for z, kmf in sorted(
        kmflist.items(), key=lambda x: int(x[0].strip("<>").split()[0])
    ):
        event_count[z] = len(df_power_storm[df_power_storm.distance_partition == z])
        kmf.plot_survival_function(ax=a, label=f"{z} km ({event_count[z]} events)")
    plt.xlim(xmin=0, xmax=36)
    plt.xlabel("Hours")
    plt.ylabel("Fraction of Outages Unrepaired")
    plt.tight_layout()
    plt.savefig("imelda_survival.png")
    plt.savefig("imelda_survival.pdf")


def plot_imelda_distance_duration(df_power_storm):
    f = plt.figure(figsize=(fig_w, fig_h))
    x, y = df_power_storm[df_power_storm.closest_flooded_road_dist < 10000][
        ["closest_flooded_road_dist", "durationh"]
    ].values.T
    plt.plot(x, y, ".", alpha=0.3)
    df_power_storm["kmcuts"] = pd.cut(
        df_power_storm.closest_flooded_road_dist, bins=np.arange(0, 22, 0.5)
    )
    plt.xlabel("Distance from Nearest Flooded Road")
    plt.ylabel("Outage Duration")
    plt.axis([0, 20, 0, 60])
    plt.tight_layout()
    plt.legend()
    plt.savefig("imelda_distance_duration.png")
    plt.savefig("imelda_distance_duration.pdf")


def plot_imelda_elastic_near_far(df_power_storm):
    
    def get_near_far_data_to_plot(df_power_storm, d):
        df = df_power_storm[
        df_power_storm.util.isin(["centerpoint_tx", "entergy_latxak"])
    ].set_index(["util", "index_copy"])
        dfclose = df[(df.closest_flooded_road_dist < d)]
        dfnotclose = df[(df.closest_flooded_road_dist >= d)]
        dfstartclose = helpers.get_outages_per_hour(dfclose, "2H")
        dfstartnotclose = helpers.get_outages_per_hour(dfnotclose, "2H")
        xc, yc = dfstartclose[["events", "ends"]].values.T
        xnc, ync = dfstartnotclose[["events", "ends"]].values.T
        res = sm.OLS(ync, sm.add_constant(xnc)).fit()
        xf = np.linspace(min(min(xc), min(xnc)), max(max(xc), max(xnc)), 1000)
        yf_nc = sm.add_constant(xf) @ res.params
        return xc, yc, xnc, ync, res, xf, yf_nc

    def plot_near_far(d, xc, yc, xnc, ync, xf, yf_nc):
        plt.plot(xf, yf_nc, color="darkblue", alpha=0.8)
        plt.plot(xnc, ync, ".", color="limegreen", label=f"flooded roads >{d} mi")
        plt.plot(xc, yc, ".", color="darkred", label=f"flooded roads <{d} mi")
        plt.legend(loc="best")

    def plot_far_only(xnc, ync, res, xf, yf_nc):
        plt.plot(xf, yf_nc, color="darkblue", alpha=0.8)
        # plt.plot(xnc,ync,'.',color="limegreen",label="flooded roads >5 mi")
        ync_pred = sm.add_constant(xnc) @ res.params
        plt.scatter(
            xnc,
            ync,
            c=ync / ync_pred,
            norm=helpers.get_divnorm(),
            cmap="coolwarm_r",
            alpha=0.85,
        )

    def plot_near_only(xc, yc, res, xf, yf_nc):
        plt.plot(xf, yf_nc, color="darkblue", alpha=0.8)
        # plt.plot(xnc,ync,'.',color="limegreen",label="flooded roads >5 mi")
        # plt.plot(xc,yc,'.',color="darkred",label="flooded roads <5 mi")
        yc_pred = sm.add_constant(xc) @ res.params
        plt.scatter(
            xc,
            yc,
            c=yc / yc_pred,
            norm=helpers.get_divnorm(),
            cmap="coolwarm_r",
            alpha=0.85,
        )

    for d in [5, 10, 20, 30]:
        xc, yc, xnc, ync, res, xf, yf_nc = get_near_far_data_to_plot(df_power_storm, d)
        for figstyle in [
            f"imelda_stress_strain_near_far_{d}mi",
            f"imelda_stress_strain_far_{d}mi",
            f"imelda_stress_strain_near_{d}mi",
        ]:
            plt.figure(figsize=(5, 4.5))
            if figstyle == f"imelda_stress_strain_near_far_{d}mi":
                plot_near_far(d, xc, yc, xnc, ync, xf, yf_nc)
            elif figstyle == f"imelda_stress_strain_far_{d}mi":
                plot_far_only(xnc, ync, res, xf, yf_nc)
            elif figstyle == f"imelda_stress_strain_near_{d}mi":
                plot_near_only(xc, yc, res, xf, yf_nc)

            plt.xlabel("Outages/2H")
            plt.ylabel("Repairs/2H")
            plt.tight_layout()
            plt.axis([0, max(max(xc), max(xnc)), 0, max(max(yc), max(ync))])
            plt.savefig(figstyle + ".png")
            plt.savefig(figstyle + ".pdf")

if __name__ == "__main__":

    df_roads = helpers.load_tx_511_imelda_df()
    df_power = helpers.load_imelda_power_df()

    imelda_start_dt = datetime.datetime(2019, 9, 17)
    imelda_end_dt = datetime.datetime(2019, 9, 27)
    df_power_storm = df_power[
        (df_power.start_timestamp > imelda_start_dt.timestamp())
        & (df_power.start_timestamp < imelda_end_dt.timestamp())
    ].reset_index()

    df_prev60 = df_power[
        (df_power.start_timestamp > datetime.datetime(2019, 7, 17).timestamp())
        & (df_power.start_timestamp < datetime.datetime(2019, 9, 17).timestamp())
    ]

    df_power_storm = helpers.add_closest_flooded_roads(df_roads, df_power_storm)
    df_roads = helpers.add_outage_details_to_df_roads(df_roads, df_power_storm)

    plot_simplified_survival_curve(df_power_storm, df_prev60)
    plot_imelda_survival(df_power_storm)
    plot_imelda_distance_duration(df_power_storm)
    plot_imelda_elastic_near_far(df_power_storm)