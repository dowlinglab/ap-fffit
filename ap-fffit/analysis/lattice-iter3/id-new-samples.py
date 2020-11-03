import sys
import gpflow
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn

from fffit.utils import (
    shuffle_and_split,
    values_scaled_to_real,
)


from fffit.models import run_gpflow_scipy
from fffit.pareto import find_pareto_set, is_pareto_efficient

sys.path.append("../")

from utils.ap import AP
from utils.prepare_samples import prepare_df_lattice
from utils.id_new_samples import rank_samples

############################# QUANTITIES TO EDIT #############################
##############################################################################

iternum = 2
gp_shuffle_seed = 9554774

##############################################################################
##############################################################################

distance_seed = 10

csv_path = "/scratch365/rdefever/ap-fffit/ap-fffit/analysis/csv/"
in_csv_names = [
    "lattice-iter" + str(i) + "-results.csv" for i in range(1, iternum+1)
]
out_csv_name = "lattice-iter" + str(iternum+1) + "-params.csv"

# Read files
df_csvs = [
    pd.read_csv(csv_path + in_csv_name, index_col=0)
    for in_csv_name in in_csv_names
]
df_csv = pd.concat(df_csvs)
df_all = prepare_df_lattice(df_csv, AP)

def main():

    ### Fit GP models
    # Create training/test set
    param_names = list(AP.param_names) + ["temperature"]
    property_names = ["sim_lattice_a", "sim_lattice_b", "sim_lattice_c"]

    models = {}
    for property_name in property_names:
        # Get train/test
        x_train, y_train, x_test, y_test = shuffle_and_split(
            df_all, param_names, property_name, shuffle_seed=gp_shuffle_seed
        )

        # Fit model
        models[property_name] = run_gpflow_scipy(
            x_train,
            y_train,
            gpflow.kernels.RBF(lengthscales=np.ones(AP.n_params + 1)),
        )

    ### Step 3: Find new parameters for simulations

    # Load large hypercube
    latin_hypercube = np.loadtxt("LHS_1e6x8.csv", delimiter=",")

    # Calculate properties
    predicted_mses = {}
    for property_name, model in models.items():
        predicted_mses[property_name] = rank_samples(
            latin_hypercube, model, AP, property_name
        )

    # Merge into single DF
    mses = predicted_mses["sim_lattice_a"].merge(
        predicted_mses["sim_lattice_b"], on=AP.param_names
    )
    mses = mses.rename(
        {"mse_x": "mse_lattice_a", "mse_y": "mse_lattice_b"}, axis="columns"
    )
    mses = mses.merge(predicted_mses["sim_lattice_c"], on=AP.param_names)
    mses = mses.rename(
        {"mse": "mse_lattice_c"}, axis="columns"
    )

    # 0.01 angstrom^2 -> 0.1 angstrom -> ~1-2%
    for mse_threshold in [0.03, 0.02, 0.01]:
        # Filter for points with all mses below threshold
        query = (f"mse_lattice_a < {mse_threshold} & "
                 f"mse_lattice_b < {mse_threshold} & "
                 f"mse_lattice_c < {mse_threshold}"
                )
        filtered = mses.query(query)

        # Plot filtered points vs. params
        g = seaborn.pairplot(filtered, vars=list(AP.param_names))
        g.set(xlim=(-0.1, 1.1), ylim=(-0.1, 1.1))
        g.savefig("figs/AP-filtered-points.pdf")

        print(f"A total of {len(filtered)} points were found with all mse < {mse_threshold}.")

    # Final mse threshold = 0.01
    # Keep the best points for each lattice
    # Sort remaining by mean MSE. Greedy search for well-spaced points.
    filtered["mean_mse"] = filtered.filter(regex="mse").mean(axis=1)
    new_points = filtered.sort_values("mse_lattice_a").iloc[[0]]
    new_points = new_points.append(filtered.sort_values("mse_lattice_b").iloc[[0]])
    new_points = new_points.append(filtered.sort_values("mse_lattice_c").iloc[[0]])
    filtered.drop(index=new_points.index, inplace=True)
    filtered = filtered.sort_values("mean_mse")

    distance = 0.7977
    discarded_points = pd.DataFrame(columns=filtered.columns)
    np.random.seed(distance_seed)

    while len(filtered > 0):
        # Grab the current best point
        new_points = new_points.append(filtered.iloc[[0]])
        # Remove anything within distance
        l1_norm = np.sum(
            np.abs(
                filtered[list(AP.param_names)].values
                - new_points[list(AP.param_names)].iloc[[-1]].values
            ),
            axis=1,
        )
        points_to_remove = np.where(l1_norm < distance)[0]
        discarded_points = discarded_points.append(
            filtered.iloc[points_to_remove]
        )
        filtered.drop(
            index=filtered.index[points_to_remove], inplace=True
        )
    print(
        f"After removing similar points, we are left with {len(new_points)} points."
    )


    # Plot new points
    g = seaborn.pairplot(new_points, vars=list(AP.param_names))
    g.set(xlim=(-0.1, 1.1), ylim=(-0.1, 1.1))
    g.savefig("figs/AP-new-points.pdf")

    # Save to CSV
    new_points.drop(
        columns=[
            "mse_lattice_a",
            "mse_lattice_b",
            "mse_lattice_c",
            "mean_mse",
        ],
        inplace=True,
    )
    new_points.to_csv(csv_path + out_csv_name)

if __name__ == "__main__":
    main()
