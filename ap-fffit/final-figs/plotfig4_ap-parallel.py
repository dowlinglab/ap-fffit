import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
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
    df = pd.read_csv("../analysis/csv/uc-lattice-final-params.csv", index_col=0)
    df = df[df.temperature == 10]
    dff = pd.read_csv("../analysis/csv/ap-final-2.csv", index_col=0)
    seaborn.set_palette('bright', n_colors=len(df))
    #data = values_real_to_scaled(df[list(AP.param_names)].values, AP.param_bounds)
    #data_f = values_real_to_scaled(dff[list(AP.param_names)].values, AP.param_bounds)
    data = df[list(AP.param_names)].values
    data_f = dff[list(AP.param_names)].values
    result_bounds = np.array([[0, 0.5], [0, 5]])
    results = values_real_to_scaled(df[["uc_mean_distance", "lattice_ape"]].values, result_bounds)
    results_f = values_real_to_scaled(dff[["uc_mean_distance", "lattice_ape"]].values, result_bounds)
    param_bounds = AP.param_bounds
    param_bounds[:4] = param_bounds[:4] #* NM_TO_ANGSTROM
    param_bounds[4:] = param_bounds[4:] #* KJMOL_TO_K

    data = np.hstack((data, results))
    data_f = np.hstack((data_f, results_f))
    bounds = np.vstack((param_bounds, result_bounds))
    print(data.shape)

    col_names = [
        r"$\sigma_{Cl}$",
        r"$\sigma_H$",
        r"$\sigma_N$",
        r"$\sigma_O$",
        r"$\epsilon_{Cl}$",
        r"$\epsilon_H$",
        r"$\epsilon_N$",
        r"$\epsilon_O$",
        r"UCMD",
        "Lattice\nMAPE",
    ]
    n_axis = len(col_names)
    assert data.shape[1] == n_axis
    x_vals = [i for i in range(n_axis)]
    
    # Create (N-1) subplots along x axis
    fig, axes = plt.subplots(1, n_axis-1, sharey=False, figsize=(12,5))
    
    # Plot each row
    for i, ax in enumerate(axes):
        for line in data:
            ax.plot(x_vals, line, alpha=0.35)
        ax.set_xlim([x_vals[i], x_vals[i+1]])
        for line in data_f:
            ax.plot(x_vals, line, alpha=1.0, linewidth=3)

    for dim, ax in enumerate(axes):
        ax.xaxis.set_major_locator(ticker.FixedLocator([dim]))
        set_ticks_for_axis(ax, bounds[dim], nticks=6)
        ax.set_xticklabels([col_names[dim]], fontsize=30)
        ax.tick_params(axis="x", pad=10)
        ax.set_ylim(-0.05,1.05)
        # Add white background behind labels
        for label in ax.get_yticklabels():
            label.set_bbox(
                dict(
                    facecolor='white',
                    edgecolor='none',
                    alpha=0.45,
                    boxstyle=mpatch.BoxStyle("round4")
                )
            )
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_linewidth(2.0)

    ax = axes[-1]
    ax.xaxis.set_major_locator(ticker.FixedLocator([n_axis-2, n_axis-1]))
    ax.set_xticklabels([col_names[-2], col_names[-1]], fontsize=24)
    ax.tick_params(axis="x", pad=14)

    # Add class II data
    #ax.plot(x_vals[-2], 0.3485/bounds[-2][1], markersize=15, color="red", marker="*", clip_on=False, zorder=200)
    
    ax = plt.twinx(axes[-1])
    ax.set_ylim(-0.05, 1.05)
    set_ticks_for_axis(ax, bounds[-1], nticks=6)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_linewidth(2.0)

    # Add class I data
    ax.plot(x_vals[-1], 1.42/bounds[-1][1], markersize=15, color="tab:purple", marker="o", clip_on=False, zorder=200)
    ax.plot(x_vals[-2], 0.156/bounds[-2][1], markersize=15, color="tab:purple", marker="o", clip_on=False, zorder=200)

    # Add class II data
    ax.plot(x_vals[-1], 3.55/bounds[-1][1], markersize=20, color="tab:red", marker="*", clip_on=False, zorder=200)
    ax.plot(x_vals[-2], 0.3485/bounds[-2][1], markersize=20, color="tab:red", marker="*", clip_on=False, zorder=200)


    # Remove space between subplots
    plt.subplots_adjust(wspace=0, bottom=0.2, left=0.05, right=0.95)
    
    fig.savefig("pdfs/fig4-ap-parallel.pdf")


def set_ticks_for_axis(ax, bounds, nticks):
    """Set the tick positions and labels on y axis for each plot

    Tick positions based on normalised data
    Tick labels are based on original data
    """
    min_val, max_val = bounds
    step = (max_val - min_val) / float(nticks-1)
    tick_labels = [round(min_val + step * i, 3) for i in range(nticks)]
    ticks = np.linspace(0, 1.0, nticks)
    ax.yaxis.set_ticks(ticks)
    ax.set_yticklabels(tick_labels, fontsize=18)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.tick_params("y", direction="inout", which="both", length=6)
    ax.tick_params("y", which="major", length=12)
    ax.tick_params("x", pad=15) 

if __name__ == "__main__":
    main()

