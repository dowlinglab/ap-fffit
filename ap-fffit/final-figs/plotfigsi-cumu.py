import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import seaborn

sys.path.append("../analysis")
from utils.ap import APConstants
AP = APConstants()

matplotlib.rc("font", family="sans-serif")
matplotlib.rc("font", serif="Arial")

############################# QUANTITIES TO EDIT #############################
##############################################################################

iternum = 4

##############################################################################
##############################################################################

csv_path = "../analysis/csv/"
in_csv_names = [
    "uc-lattice-iter" + str(i) + "-results.csv" for i in range(1, iternum + 1)
]

# Read files
dfs_csv = [
    pd.read_csv(csv_path + in_csv_name, index_col=0)
    for in_csv_name in in_csv_names
]

dfs = []
for df in dfs_csv:
    means = df.groupby(by=list(AP.param_names))["lattice_ape"].mean()
    df = df.join(means, on=list(AP.param_names), rsuffix="_mean")
    df = df[df.temperature == 10]
    dfs.append(df)

def main():

    # Plot MAPE sorted by each property
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))
    name = "uc_mean_distance"
    for iiter in range(iternum):
        ax1.plot(
            dfs[iiter].sort_values(name)[name],
            np.arange(1, 251,1),
            '-o',
            markersize=6,
            alpha=0.6,
            label=f"Iteration {iiter+1}",
        )


    ax1.xaxis.set_major_locator(MultipleLocator(1))
    ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.tick_params("both", direction="in", which="both", length=2, labelsize=14)
    ax1.set_xlabel(r"UCMD ($\mathrm{\AA}$)", fontsize=18, labelpad=10)
    ax1.set_ylabel(r"$N_\mathrm{cumu.}$ param. sets", fontsize=18, labelpad=14)
    ax1.tick_params("both", which="major", length=4)
    ax1.xaxis.set_ticks_position("both")
    ax1.yaxis.set_ticks_position("both")
    ax1.text(0.8, 0.8, "(a)", fontsize=22, transform=ax1.transAxes)

    ax1ins = inset_axes(ax1, width="100%", height="100%",
            bbox_to_anchor=(0.4, 0.10, 0.45, 0.45),
            bbox_transform=ax1.transAxes, loc=3)
    for iiter in range(iternum):
        ax1ins.plot(
            dfs[iiter].sort_values(name)[name],
            np.arange(1, 251,1),
            '-o',
            markersize=3,
            alpha=0.6,
            label=f"Iteration {iiter+1}",
        )

    ax1ins.set_xlim(-0.05, 0.35)
    ax1ins.set_ylim(-5, 105)
    ax1ins.xaxis.set_major_locator(MultipleLocator(0.15))
    ax1ins.tick_params("both", direction="in", which="both")

    name = "lattice_ape_mean"
    for iiter in range(iternum):
        ax2.plot(
            dfs[iiter].sort_values(name)[name],
            np.arange(1, 251,1),
            '-o',
            markersize=6,
            alpha=0.6,
            label=f"Iteration {iiter+1}",
        )

    ax2.xaxis.set_major_locator(MultipleLocator(5))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.tick_params("both", direction="in", which="both", length=2, labelsize=14)
    ax2.set_xlabel(r"Lattice MAPE", fontsize=18, labelpad=10)
    ax2.set_ylabel(r"$N_\mathrm{cumu.}$ param. sets", fontsize=18, labelpad=14)
    ax2.tick_params("both", which="major", length=4)
    ax2.xaxis.set_ticks_position("both")
    ax2.yaxis.set_ticks_position("both")
    ax2.text(0.8, 0.8, "(b)", fontsize=22, transform=ax2.transAxes)

    ax2ins = inset_axes(ax2, width="100%", height="100%",
            bbox_to_anchor=(0.52, 0.08, 0.4, 0.4),
            bbox_transform=ax2.transAxes, loc=3)
    for iiter in range(iternum):
        ax2ins.plot(
            dfs[iiter].sort_values(name)[name],
            np.arange(1, 251,1),
            '-o',
            markersize=3,
            alpha=0.6,
            label=f"Iteration {iiter+1}",
        )

    ax2ins.set_xlim(-0.1, 2.1)
    ax2ins.set_ylim(-5, 105)
    ax2ins.xaxis.set_major_locator(MultipleLocator(1))
    ax2ins.tick_params("both", direction="in", which="both")

    plt.subplots_adjust(left = 0.15, right=0.95, top=0.8, bottom=0.2, wspace=0.5)
    ax1.legend(fontsize=16, loc=(-0.4,1.07), ncol=4, columnspacing=1, handletextpad=0.5)
    fig.savefig("pdfs/figsi-ap-cumu.pdf")

if __name__ == "__main__":
    main()
