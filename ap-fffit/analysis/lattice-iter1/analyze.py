import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn

from fffit.utils import (
    values_real_to_scaled,
    values_scaled_to_real,
)

sys.path.append("../")

from utils.ap import AP
from utils.prepare_samples import prepare_df_lattice
from utils.prepare_samples import prepare_df_lattice_errors
from utils.plot import plot_property, render_mpl_table

import matplotlib._color_data as mcd

############################# QUANTITIES TO EDIT #############################
##############################################################################

iternum = 1

##############################################################################
##############################################################################

csv_path = "/scratch365/rdefever/ap-fffit/ap-fffit/analysis/csv/"
in_csv_names = [
    "lattice-iter" + str(i) + "-results.csv" for i in range(1, iternum + 1)
]

# Read files
df_csvs = [
    pd.read_csv(csv_path + in_csv_name, index_col=0)
    for in_csv_name in in_csv_names
]
df_csv = pd.concat(df_csvs)
df_all = prepare_df_lattice(df_csv, AP)


def main():

    # Basic plots to view output
    plot_property(
        df_all,
        "lattice_a",
        AP.lattice_a_bounds,
        axis_name="a [angstrom]",
    )
    plot_property(
        df_all,
        "lattice_b",
        AP.lattice_b_bounds,
        axis_name="b [angstrom]",
    )
    plot_property(
        df_all,
        "lattice_c",
        AP.lattice_c_bounds,
        axis_name="c [angstrom]",
    )

    # Create a dataframe with one row per parameter set
    df_paramsets = prepare_df_lattice_errors(df_all, AP)

    # Plot MSE for each property
    fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(5, 10))
    fig.suptitle("Mean square errors", y=1.05, x=0.6)
    ax1.plot(
        range(df_paramsets.shape[0]),
        df_paramsets["mse_lattice_a"],
        label="a",
        color="black",
        marker="o",
        linestyle="",
        ms=3,
    )
    ax1.set_title("Lattice a")
    ax1.set_ylabel(r"angstrom$^2$")
    ax2.plot(
        range(df_paramsets.shape[0]),
        df_paramsets["mse_lattice_b"],
        label="b",
        color="black",
        marker="o",
        linestyle="",
        ms=3,
    )
    ax2.set_title("Lattice b")
    ax2.set_ylabel(r"angstrom$^2$")
    ax3.plot(
        range(df_paramsets.shape[0]),
        df_paramsets["mse_lattice_c"],
        label="c",
        color="black",
        marker="o",
        linestyle="",
        ms=3,
    )
    ax3.set_title("Lattice c")
    ax3.set_ylabel(r"angstrom$^2$")
    ax3.set_xlabel("Index (Unsorted)")
    fig.tight_layout()
    fig.savefig("figs/mse.png", dpi=300)

    # Plot MAPE with all properties
    fig, ax = plt.subplots()
    names = {
        "mape_lattice_a": "Lattice a",
        "mape_lattice_b": "Lattice b",
        "mape_lattice_c": "Lattice c",
    }
    for name, label in names.items():
        ax.plot(
        range(df_paramsets.shape[0]), df_paramsets[name], label=label, marker="o", linestyle="", ms=3
        )
    ax.set_xlabel("Index (Unsorted)", fontsize=16, labelpad=15)
    ax.set_ylabel("Mean abs. % error", fontsize=16, labelpad=15)
    ax.tick_params(axis="both", labelsize=12)
    fig.legend(fontsize=12)
    fig.tight_layout()
    fig.savefig("figs/mape_unsorted.png", dpi=300)

    # Plot MAPE sorted by average
    fig, ax = plt.subplots()
    df_paramsets["mape_avg"] = df_paramsets.filter(regex=("mape")).mean(axis=1)
    names = {
        "mape_lattice_a": "Lattice a",
        "mape_lattice_b": "Lattice b",
        "mape_lattice_c": "Lattice c",
    }
    for name, label in names.items():
        ax.plot(
            range(df_paramsets.shape[0]),
            df_paramsets.sort_values("mape_avg")[name],
            label=label,
        )
    ax.set_xlabel("Index Sorted by Mean MAPE", fontsize=16, labelpad=15)
    ax.set_ylabel("Mean abs. % error", fontsize=16, labelpad=15)
    ax.tick_params(axis="both", labelsize=12)
    fig.legend(fontsize=12)
    fig.tight_layout()
    fig.savefig("figs/mape_sorted.png", dpi=300)


    # Generate colors for param sets
    phi = np.linspace(0, 2 * np.pi, len(df_paramsets))
    rgb_cycle = np.vstack(
        (  # Three sinusoids
            0.5 * (1.0 + np.cos(phi)),  # scaled to [0,1]
            0.5 * (1.0 + np.cos(phi + 2 * np.pi / 3)),  # 120° phase shifted.
            0.5 * (1.0 + np.cos(phi - 2 * np.pi / 3)),
        )
    ).T  # Shape = (df_paramsets,3)

    # Plot Lattice a/b/c
    fig, ax = plt.subplots()
    ax.set_xlabel("Temperature [K]", fontsize=16, labelpad=15)
    ax.set_ylabel("Lattice constant [angstrom]", fontsize=16, labelpad=15)
    ax.tick_params(axis="both", labelsize=12)
    ax.set_xlim(-10, 320)

    for temp in AP.expt_lattice_a.keys():
        ax.scatter(
            np.tile(temp, len(df_paramsets)),
            df_paramsets.filter(regex=(f"lattice_a_{float(temp):.0f}K")),
            c=rgb_cycle,
            s=60,
            alpha=0.5,
        )
    ax.scatter(
        AP.expt_lattice_a.keys(),
        AP.expt_lattice_a.values(),
        color="black",
        marker="x",
        s=80,
    )

    fig.tight_layout()
    fig.savefig("figs/lattice_a.png", dpi=300)

    fig, ax = plt.subplots()
    ax.set_xlabel("Temperature [K]", fontsize=16, labelpad=15)
    ax.set_ylabel("Lattice constant [angstrom]", fontsize=16, labelpad=15)
    ax.tick_params(axis="both", labelsize=12)
    ax.set_xlim(-10, 320)

    for temp in AP.expt_lattice_b.keys():
        ax.scatter(
            np.tile(temp, len(df_paramsets)),
            df_paramsets.filter(regex=(f"lattice_b_{float(temp):.0f}K")),
            c=rgb_cycle,
            s=60,
            alpha=0.5,
        )
    ax.scatter(
        AP.expt_lattice_b.keys(),
        AP.expt_lattice_b.values(),
        color="black",
        marker="x",
        s=80,
    )

    fig.tight_layout()
    fig.savefig("figs/lattice_b.png", dpi=300)

    fig, ax = plt.subplots()
    ax.set_xlabel("Temperature [K]", fontsize=16, labelpad=15)
    ax.set_ylabel("Lattice constant [angstrom]", fontsize=16, labelpad=15)
    ax.tick_params(axis="both", labelsize=12)
    ax.set_xlim(-10, 320)

    for temp in AP.expt_lattice_c.keys():
        ax.scatter(
            np.tile(temp, len(df_paramsets)),
            df_paramsets.filter(regex=(f"lattice_c_{float(temp):.0f}K")),
            c=rgb_cycle,
            s=60,
            alpha=0.5,
        )
    ax.scatter(
        AP.expt_lattice_c.keys(),
        AP.expt_lattice_c.values(),
        color="black",
        marker="x",
        s=80,
    )

    fig.tight_layout()
    fig.savefig("figs/lattice_c.png", dpi=300)


    # Save tables of data sorted by sum of squares results
    render_mpl_table(
        df_paramsets.sort_values("mape_lattice_a").head(),
        out_name="table_sort_lattice_a",
    )
    render_mpl_table(
        df_paramsets.sort_values("mape_lattice_b").head(),
        out_name="table_sort_lattice_b",
    )
    render_mpl_table(
        df_paramsets.sort_values("mape_lattice_c").head(),
        out_name="table_sort_lattice_c",
    )

    # Pair plots to look at correlation of errors
    grid = seaborn.PairGrid(
        data=df_paramsets,
        vars=[
            "mape_lattice_a",
            "mape_lattice_b",
            "mape_lattice_c",
        ],
    )
    grid = grid.map_lower(plt.scatter)
    grid = grid.map_diag(plt.hist, bins=10, edgecolor="k")
    for i, j in zip(*np.triu_indices_from(grid.axes, 1)):
        grid.axes[i, j].set_visible(False)

    grid.savefig("figs/mape_pairs.png")


if __name__ == "__main__":
    main()
