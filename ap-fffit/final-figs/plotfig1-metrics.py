import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import seaborn


matplotlib.rc("font", family="sans-serif")
matplotlib.rc("font", serif="Arial")

df_dom = pd.read_csv("dat/fig1_DominatedPoints.csv", names=["metric1", "metric2"], sep="\s+")
df_dom = df_dom.sample(n=10000)
df_nodom = pd.read_csv("dat/fig1_NonDominatedPoints.csv", names=["metric1", "metric2"], sep="\s+")

def main():

    seaborn.set_palette("Paired")

    # Plot
    fig, ax = plt.subplots()

    ax.scatter(
        df_dom["metric1"],
        df_dom["metric2"],
        label="Dominated",
        alpha=0.08,
        s=90,
        c="C1",
    )
    ax.scatter(
        df_nodom["metric1"],
        df_nodom["metric2"],
        marker="+",
        label="Non-dominated",
        alpha=0.6,
        s=130,
        c="C5",
    )

    ax.tick_params("both", direction="in", which="both", length=4, labelsize=26, pad=10)
    ax.tick_params("both", which="major", length=8)
    ax.xaxis.set_ticks_position("both")
    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-2, 22)

    ax.set_xlabel(r"Metric 1", fontsize=28, labelpad=20)
    ax.set_ylabel(r"Metric 2", fontsize=28, labelpad=10)
    leg = ax.legend(
        fontsize=24,
        handletextpad=-0.5,
        borderaxespad=0.2,
        borderpad=0.2,
        loc="upper left",
        bbox_to_anchor=(0.01, 0.95)
    )
    for lh in leg.legendHandles:
            lh.set_alpha(0.7)

    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable="box")
    fig.tight_layout()
    fig.savefig("pdfs/fig1-metrics.pdf")
    fig.savefig("pdfs/fig1-metrics.png")


if __name__ == "__main__":
    main()
