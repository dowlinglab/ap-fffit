import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import seaborn

sys.path.append("../analysis/")
from utils.ap import APConstants

AP = APConstants()

matplotlib.rc("font", family="sans-serif")
matplotlib.rc("font", serif="Arial")

df_dom = pd.read_csv("../analysis/csv/uc-lattice-final-params.csv", index_col=0)
df_nondom = pd.read_csv("../analysis/csv/ap-final-2.csv", index_col=0)

means = df_dom.groupby(by=list(AP.param_names))["lattice_ape"].mean()
df_dom = df_dom.join(means, on=list(AP.param_names), rsuffix="_mean")
df_dom = df_dom[df_dom.temperature == 10]

def main():

    seaborn.set_palette("Paired")

    # Plot
    fig, ax = plt.subplots()

    ax.scatter(
        df_dom["uc_mean_distance"],
        df_dom["lattice_ape_mean"],
        label="Dominated",
        alpha=0.6,
        s=90,
        c="C1",
    )
    ax.scatter(
        df_nondom["uc_mean_distance"],
        df_nondom["Lattice_MAPE"],
        label="Non-dominated",
        alpha=1.0,
        s=90,
        c="C5",
    )

    ax.tick_params("both", direction="in", which="both", length=4, labelsize=18, pad=6)
    ax.tick_params("both", which="major", length=8)
    ax.xaxis.set_ticks_position("both")
    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.set_xlim(-0.02, 0.18)
    ax.set_ylim(-0.1, 1.6)

    ax.set_xlabel(r"UCMD ($\mathrm{\AA}$)", fontsize=22, labelpad=12)
    ax.set_ylabel(r"Lattice MAPE", fontsize=22, labelpad=12)
    leg = ax.legend(
        fontsize=16,
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
    fig.savefig("pdfs/figsi-pareto.pdf")


if __name__ == "__main__":
    main()
