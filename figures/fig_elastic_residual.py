"""Plot elastic residual figure.

Figure shows deviations from number of repairs executed in given time window
compared to expected number if linear ("elastic") behavior had continued.
"""
# pylint:disable=invalid-name
from sys import argv
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
from pandas.plotting import register_matplotlib_converters
import helpers

sb.set("paper", font_scale=1.9)
register_matplotlib_converters()


fig_w, fig_h = (8, 8)


def make_elastic_residual_plot_semilogx_deviation_colored(
    dfres: pd.DataFrame,
) -> None:
    """Make plot differences between elastic repairs and actual repairs.

    The plot will be on a semilogx scale and the 2H periods colored
    by the extent of their deviation from the elastic behavior.

    Args:
        dfres (pd.DataFrame): residuals dataframe
    """
    # make semilogx residual plot
    plt.figure(figsize=(fig_w, fig_h))
    sb.scatterplot(
        x="X",
        y="residual",
        hue="residual_fraction",
        hue_norm=(0, 2),
        data=dfres.reset_index().rename(columns=dict(level_0="Utility")),
        palette="coolwarm_r",
        legend=False,
        linewidth=0,
        alpha=0.9,
    )
    bin_x, mean, std = helpers.get_weighted_errorbar_line(dfres.X, dfres.residual)
    plt.errorbar(bin_x, mean, yerr=std, color="gray", lw=1.5)
    plt.xscale("log")
    plt.xlabel("Outage Count")
    # plt.ylabel("actual repairs - elastic repairs")
    plt.ylabel("Elastic Residual (observed repairs - predicted repairs)/2H")
    # plt.legend(ncol=int(round(np.sqrt(until_util)))-1,loc="best")
    plt.legend(ncol=2, loc="best")
    plt.tight_layout()
    plt.savefig(f"elastic_residual_top_{until_util}_rolling.png")
    plt.savefig(f"elastic_residual_top_{until_util}_rolling.pdf")


def make_elastic_residual_plot_loglog_util_colored(dfres: pd.DataFrame) -> None:
    """Make a residual fraction plot.

    Residual fraction is the ratio of the actual repairs to elastic repairs.

    Args:
        dfres (pd.DataFrame): residuals data frame
    """
    # make loglog residual fraction plot
    plt.figure(figsize=(fig_w, fig_h))
    sb.scatterplot(
        x="X",
        y="residual_fraction",
        hue="Utility",
        data=dfres.reset_index().rename(columns=dict(level_0="Utility")),
    )

    bin_x, mean, std = helpers.get_weighted_errorbar_line(
        dfres.X, dfres.residual_fraction
    )
    plt.errorbar(bin_x, mean, yerr=std, color="gray", lw=1.5)

    plt.xscale("log")
    plt.xlabel("Outage Count")
    # plt.xlabel("outage count")
    plt.ylabel("actual repairs / elastic repairs")
    plt.yscale("log")
    plt.axis(helpers.non_zero_bounds(dfres.X, dfres.residual_fraction))
    plt.legend(ncol=2, loc="best")
    plt.tight_layout()
    plt.savefig(f"elastic_residual_loglog_top_{until_util}_rolling.png")
    plt.savefig(f"elastic_residual_loglog_top_{until_util}_rolling.pdf")


if __name__ == "__main__":

    until_util = int(argv[1]) if len(argv) > 1 else 30
    dfstart_fname = argv[2] if len(argv) > 2 else None
    rolling = True

    dfstart = helpers.load_dfstart(dfstart_fname)
    mostevents = (
        dfstart.groupby(level=0).sum().starts.sort_values(ascending=False).index
    )

    residuals_dataframe = pd.concat(
        dict(
            (util, helpers.get_Y_Yhat(dfstart, util, rolling))
            for util in mostevents[:until_util]
        )
    )
    make_elastic_residual_plot_semilogx_deviation_colored(residuals_dataframe)
    make_elastic_residual_plot_loglog_util_colored(residuals_dataframe)
