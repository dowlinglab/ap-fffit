import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import seaborn

sys.path.append("../analysis")

from fffit.utils import values_real_to_scaled
from utils.ap import APConstants
from matplotlib import ticker

AP = APConstants()

matplotlib.rc("font", family="sans-serif")
matplotlib.rc("font", serif="Arial")

NM_TO_ANGSTROM = 10
K_B = 0.008314 # J/MOL K
KJMOL_TO_K = 1.0 / K_B

def main():

    df = pd.read_csv("../analysis/csv/uc-lattice-iter4-results.csv", index_col=0)

    seaborn.set_palette(None, n_colors=6)

    fig, ax = plt.subplots()
    ax.scatter(
        range(len(df.loc[df["temperature"] == 10])),
        df.loc[df["temperature"] == 10].sort_values("uc_mean_distance")["uc_mean_distance_Cl"],
        label="UCMD Cl",
        alpha=0.4,
        s=20,
        marker="v",
    )
    ax.scatter(
        range(len(df.loc[df["temperature"] == 10])),
        df.loc[df["temperature"] == 10].sort_values("uc_mean_distance")["uc_mean_distance_O"],
        label="UCMD O",
        alpha=0.4,
        s=16,
        marker="s",
    )
    ax.scatter(
        range(len(df.loc[df["temperature"] == 10])),
        df.loc[df["temperature"] == 10].sort_values("uc_mean_distance")["uc_mean_distance_N"],
        label="UCMD N",
        alpha=0.4,
        marker="P",
        s=16,
    )
    ax.scatter(
        range(len(df.loc[df["temperature"] == 10])),
        df.loc[df["temperature"] == 10].sort_values("uc_mean_distance")["uc_mean_distance_H"],
        label="UCMD H",
        alpha=0.4,
        s=20,
        marker="*",
    )
    ax.plot(
        range(len(df.loc[df["temperature"] == 10])),
        df.loc[df["temperature"] == 10].sort_values("uc_mean_distance")["uc_mean_distance"],
        label="UCMD",
        linewidth=4,
        alpha=0.8,
        color="gray",
    )
    ax2 = ax.twinx()
    ax2.scatter(
        range(len(df.loc[df["temperature"] == 10])),
        df.loc[df["temperature"] == 10].sort_values("uc_mean_distance")["abs(HB3-HB4)"],
        label=r"abs$(\Delta_\mathrm{hbond})$",
        color="black",
        alpha=0.4,
        s=20,
    )

    ax.set_xlabel("Sorted parameter index", fontsize=16, labelpad=15)
    ax.set_ylabel(r"UCMD [$\mathrm{\AA}$]", fontsize=16, labelpad=15)
    ax2.set_ylabel(r"abs$(\Delta_\mathrm{hbond})$ $[\mathrm{\AA}]$", fontsize=16, labelpad=15)

    ax2.set_ylim(0.0,0.004)

    ax.tick_params("both", direction="in", which="both", length=4, labelsize=14, pad=5)
    ax2.tick_params("both", direction="in", which="both", length=4, labelsize=14, pad=5)
    ax.tick_params("both", which="major", length=8)
    ax2.tick_params("both", which="major", length=8)
    ax.xaxis.set_ticks_position("both")
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.set_ylim(0.0, 0.4)

    fig.subplots_adjust(left=0.2, bottom=0.2, top=0.8, right=0.8)

    fig.legend(fontsize=12, loc="lower right", bbox_to_anchor=(0.9,1), bbox_transform=ax.transAxes, ncol=2)
    fig.savefig("pdfs/fig8-ap-tradeoffs.pdf")


if __name__ == "__main__":
    main()
