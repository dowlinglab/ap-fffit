import os
import sys
import gpflow
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn

from sklearn import svm
from sklearn import metrics

from fffit.utils import (
    shuffle_and_split,
    values_scaled_to_real,
    values_real_to_scaled,
)


from fffit.models import run_gpflow_scipy
from fffit.pareto import find_pareto_set, is_pareto_efficient

sys.path.append("../")

from utils.ap import AP
from utils.prepare_samples import prepare_df

############################# QUANTITIES TO EDIT #############################
##############################################################################

iternum = 2
gp_shuffle_seed = 588654
clf_shuffle_seed = 19485
distance_seed = 15

##############################################################################
##############################################################################

temperatures = [10, 78, 298]
ucmd_clf_threshold = 0.8
ucmd_next_itr_threshold = 0.35
lattice_mape_next_itr_threshold = 2.5

csv_path = "/scratch365/bbefort/ap-fffit/ap-fffit/analysis/csv/"
in_csv_names = [
    "uc-lattice-iter" + str(i) + "-results.csv" for i in range(1, iternum+1)
]
out_csv_name = "uc-lattice-iter" + str(iternum+1) + "-params.csv"

# Read files
df_csvs = [
    pd.read_csv(csv_path + in_csv_name, index_col=0)
    for in_csv_name in in_csv_names
]
df_csv = pd.concat(df_csvs)
df_all = prepare_df(df_csv, AP)


def main():

    # Make sure we have a folder for figures
    try:
        os.mkdir("figs")
    except FileExistsError:
        pass

    ###########################################################
    ####################   Fit GP models    ###################
    ###########################################################
    ## UCMD Model
    param_names = list(AP.param_names)
    property_name = "uc_mean_distance"
    # Only train on 10 K data that meets ucmd clf threshold
    df_ucmd = df_all.loc[(df_all["temperature"] == 10) & (df_all["uc_mean_distance"] < ucmd_clf_threshold)]
    # Get train/test
    x_train, y_train, x_test, y_test = shuffle_and_split(
        df_ucmd, param_names, property_name, shuffle_seed=gp_shuffle_seed
    )
    # Fit model
    ucmd_gp = run_gpflow_scipy(
        x_train,
        y_train,
        gpflow.kernels.Matern32(lengthscales=np.ones(AP.n_params)),
    )

    ## Lattice APE Model
    param_names = list(AP.param_names) + ["scaled_temperature"]
    property_name = "lattice_ape"
    # Get train/test
    x_train, y_train, x_test, y_test = shuffle_and_split(
        df_all, param_names, property_name, shuffle_seed=gp_shuffle_seed
    )
    # Fit model
    lattice_gp = run_gpflow_scipy(
        x_train,
        y_train,
        gpflow.kernels.Matern32(lengthscales=np.ones(AP.n_params + 1)),
    )

    ###########################################################
    ##################   Train classifers    ##################
    ###########################################################

    x_train, y_train, x_test, y_test = shuffle_and_split(
        df_all.loc[df_all["temperature"] == 10],
        list(AP.param_names),
        "uc_mean_distance",
        shuffle_seed=clf_shuffle_seed,
    )

    y_train = np.array(y_train < ucmd_clf_threshold, dtype=np.int32)
    y_test = np.array(y_test < ucmd_clf_threshold, dtype=np.int32)
    ucmd_clf = svm.SVC(kernel="rbf")
    ucmd_clf.fit(x_train, y_train)
    y_pred = ucmd_clf.predict(x_train)
    print("Training accuracy:",metrics.accuracy_score(y_train, y_pred))
    y_pred = ucmd_clf.predict(x_test)
    print("Testing accuracy:",metrics.accuracy_score(y_test, y_pred))
    print(metrics.confusion_matrix(y_test, y_pred))


    ###########################################################
    ###################   Find new points   ###################
    ###########################################################

    # Load large hypercube
    latin_hypercube = np.loadtxt("LHS_1e6x8.csv", delimiter=",")

    # Apply UCMD classifier
    ucmd_pred = ucmd_clf.predict(latin_hypercube)

    # Predict UCMD with GP model
    gp_means_ucmd, gp_vars_ucmd = ucmd_gp.predict_f(latin_hypercube)
    # Predict Lattice APE with GP model at each temperature
    all_errs = np.empty(shape=(latin_hypercube.shape[0],len(temperatures)))
    col_idx = 0
    for temperature in temperatures:
        scaled_temperature = values_real_to_scaled(temperature, AP.temperature_bounds)
        xx = np.hstack((latin_hypercube, np.tile(scaled_temperature,(latin_hypercube.shape[0],1))))
        gp_means_ape, gp_vars_ape = lattice_gp.predict_f(xx)
        all_errs[:, col_idx] = gp_means_ape[:,0]
        col_idx += 1
    # Compute MAPE across all three temperatures
    mean_errs = np.mean(all_errs, axis=1)

    ## Save all the results to a dataframe
    LH_results = pd.DataFrame(latin_hypercube, columns=AP.param_names)
    LH_results["ucmd_clf"] = ucmd_pred.astype(np.bool_)
    LH_results["ucmd"] = gp_means_ucmd.numpy()
    LH_results["lattice_mape"] = mean_errs

    # Only take points where structure classifier is satisifed
    LH_results_pass_ucmd_clf = LH_results.loc[LH_results.ucmd_clf == True]
    costs = LH_results_pass_ucmd_clf[["ucmd", "lattice_mape"]].to_numpy()

    # Find pareto efficient points
    result, pareto_points, dominated_points = find_pareto_set(
        costs, is_pareto_efficient
    )
    LH_results_pass_ucmd_clf["is_pareto"] = result

    # Plot pareto points vs. costs
    g = seaborn.pairplot(
        LH_results_pass_ucmd_clf,
        vars=["ucmd", "lattice_mape"],
        hue="is_pareto",
    )
    g.savefig("figs/pareto-mses.png", dpi=300)

    # Plot pareto points vs. params
    g = seaborn.pairplot(LH_results_pass_ucmd_clf, vars=list(AP.param_names), hue="is_pareto")
    g.set(xlim=(-0.1, 1.1), ylim=(-0.1, 1.1))
    g.savefig("figs/pareto-params.png", dpi=300)

    # For next iteration: 1. All non-dominated points that meet the thresholds
    #                     2. "Separated" dominated points that meet the thresholds

    next_iteration_points = LH_results_pass_ucmd_clf.loc[
            (LH_results_pass_ucmd_clf.is_pareto == True) &
            (LH_results_pass_ucmd_clf.ucmd < ucmd_next_itr_threshold) &
            (LH_results_pass_ucmd_clf.lattice_mape < lattice_mape_next_itr_threshold)
    ]
    dominated_points = LH_results_pass_ucmd_clf.loc[
            (LH_results_pass_ucmd_clf.is_pareto == False) &
            (LH_results_pass_ucmd_clf.ucmd < ucmd_next_itr_threshold) &
            (LH_results_pass_ucmd_clf.lattice_mape < lattice_mape_next_itr_threshold)
    ]


    print(f"{len(LH_results_pass_ucmd_clf)} points meet the ucmd classifier.")
    print(f"{len(LH_results_pass_ucmd_clf[LH_results_pass_ucmd_clf.is_pareto == True])} are non-dominated.")
    print(f"{len(dominated_points)} are dominated with ucmd < {ucmd_next_itr_threshold} and lattice_mape < {lattice_mape_next_itr_threshold}")

    removal_distance = 1.034495
    np.random.seed(distance_seed)
    discarded_points = pd.DataFrame(columns=dominated_points.columns)

    while len(dominated_points > 0):
        # Shuffle the top parameter sets
        dominated_points = dominated_points.sample(frac=1)

        # Select one off the top
        # Note: here we use a random one rather than the "best" one; we don't have
        # confidence that the GP models are more accurate than our thresholds
        next_iteration_points = next_iteration_points.append(dominated_points.iloc[[0]])

        # Remove anything within given distance
        l1_norm = np.sum(
            np.abs(
                dominated_points[list(AP.param_names)].values
                - next_iteration_points[list(AP.param_names)].iloc[[-1]].values
                ),
                axis=1,
        )

        points_to_remove = np.where(l1_norm < removal_distance)[0]
        discarded_points = discarded_points.append(dominated_points.iloc[points_to_remove])
        dominated_points.drop(index=dominated_points.index[points_to_remove], inplace=True)

    print(f"After removing similar points, we are left with {len(next_iteration_points)} final top points.")

    next_iteration_points = next_iteration_points[:250]
    next_iteration_points.drop(columns=["ucmd","lattice_mape","is_pareto"], inplace=True)

    # Plot new points
    g = seaborn.pairplot(next_iteration_points, vars=list(AP.param_names))
    g.set(xlim=(-0.1, 1.1), ylim=(-0.1, 1.1))
    g.savefig("figs/new-points-params.png", dpi=300)

    # Save the final new parameters
    next_iteration_points.to_csv(csv_path + out_csv_name)


if __name__ == "__main__":
    main()


