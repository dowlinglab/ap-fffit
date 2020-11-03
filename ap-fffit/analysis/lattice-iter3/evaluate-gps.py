import sys
import gpflow
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages

from fffit.utils import (
    shuffle_and_split,
    values_real_to_scaled,
    values_scaled_to_real,
    variances_scaled_to_real,
)

from fffit.plot import (
    plot_model_performance,
    plot_slices_temperature,
    plot_slices_params,
    plot_model_vs_test,
)

from fffit.models import run_gpflow_scipy

sys.path.append("../")

from utils.ap import AP
from utils.prepare_samples import prepare_df_lattice

pdf = PdfPages('figs/gp_models_eval.pdf')
mpl.rc('figure', max_open_warning = 0)

############################# QUANTITIES TO EDIT #############################
##############################################################################

iternum = 3
gp_shuffle_seed = 59457837

##############################################################################
##############################################################################

csv_path = "/scratch365/rdefever/ap-fffit/ap-fffit/analysis/csv/"
in_csv_names = [
    "lattice-iter" + str(i) + "-results.csv" for i in range(1, iternum+1)
]

# Read files
df_csvs = [
    pd.read_csv(csv_path + in_csv_name, index_col=0)
    for in_csv_name in in_csv_names
]
df_csv = pd.concat(df_csvs)
df_all = prepare_df_lattice(df_csv, AP)


def main():

    ### Fit GP Model to lattice a
    param_names = list(AP.param_names) + ["temperature"]
    property_name = "sim_lattice_a"
    x_train, y_train, x_test, y_test = shuffle_and_split(
        df_all, param_names, property_name, shuffle_seed=gp_shuffle_seed, fraction_train=0.8
    )
    print("Fitting models a")
    # Fit model
    models = {}
    models["RBF"] = run_gpflow_scipy(
        x_train,
        y_train,
        gpflow.kernels.RBF(lengthscales=np.ones(AP.n_params + 1)),
    )
    models["Matern32"] = run_gpflow_scipy(
        x_train,
        y_train,
        gpflow.kernels.Matern32(lengthscales=np.ones(AP.n_params + 1)),
    )
    models["Matern52"] = run_gpflow_scipy(
        x_train,
        y_train,
        gpflow.kernels.Matern52(lengthscales=np.ones(AP.n_params + 1)),
    )

    # Plot model performance on train and test points
    pdf.savefig(plot_model_performance(models, x_train, y_train, AP.lattice_a_bounds))
    pdf.savefig(plot_model_performance(models, x_test, y_test, AP.lattice_a_bounds))

    print("Plotting models a...")
    # Plot temperature slices
    figs = plot_slices_temperature(
        models,
        AP.n_params,
        AP.temperature_bounds,
        AP.lattice_a_bounds,
        plot_bounds=[-10.0,320.0],
        property_name="Lattice a [angstrom]",
    )

    for fig in figs:
        pdf.savefig(fig)
    del figs

    # Plot parameter slices
    for param_name in AP.param_names:
        figs = plot_slices_params(
            models,
            param_name,
            AP.param_names,
            298,
            AP.temperature_bounds,
            AP.lattice_a_bounds,
            property_name="Lattice a [angstrom]",
        )
        for fig in figs:
            pdf.savefig(fig)
        del figs

    # Loop over test params
    for test_params in x_test[:,:AP.n_params]:
        train_points = []
        test_points = []
        # Locate rows where parameter set == test parameter set
        matches = np.unique(np.where((df_all[list(AP.param_names)] == test_params).all(axis=1))[0])
        # Loop over all matches -- these will be different temperatures
        for match in matches:
            # If the match (including T) is in the test set, then append to test points
            if np.where((df_all.values[match,:AP.n_params+1] == x_test[:,:AP.n_params+1]).all(axis=1))[0].shape[0] == 1:
                test_points.append([df_all["temperature"].iloc[match],df_all[property_name].iloc[match]])
            # Else append to train points
            else:
                train_points.append([df_all["temperature"].iloc[match],df_all[property_name].iloc[match]])

        pdf.savefig(
            plot_model_vs_test(
                models,
                test_params,
                np.asarray(train_points),
                np.asarray(test_points),
                AP.temperature_bounds,
                AP.lattice_a_bounds,
                plot_bounds=[-10.0, 320.0],
                property_name="Lattice a [angstrom]"
            )
        )

    ### Fit GP Model to lattice b
    param_names = list(AP.param_names) + ["temperature"]
    property_name = "sim_lattice_b"
    x_train, y_train, x_test, y_test = shuffle_and_split(
        df_all, param_names, property_name, shuffle_seed=gp_shuffle_seed, fraction_train=0.8
    )

    # Fit model
    print("Fitting models b")
    models = {}
    models["RBF"] = run_gpflow_scipy(
        x_train,
        y_train,
        gpflow.kernels.RBF(lengthscales=np.ones(AP.n_params + 1)),
    )
    models["Matern32"] = run_gpflow_scipy(
        x_train,
        y_train,
        gpflow.kernels.Matern32(lengthscales=np.ones(AP.n_params + 1)),
    )
    models["Matern52"] = run_gpflow_scipy(
        x_train,
        y_train,
        gpflow.kernels.Matern52(lengthscales=np.ones(AP.n_params + 1)),
    )

    # Plot model performance on train and test points
    pdf.savefig(plot_model_performance(models, x_train, y_train, AP.lattice_b_bounds))
    pdf.savefig(plot_model_performance(models, x_test, y_test, AP.lattice_b_bounds))

    # Plot temperature slices
    print("Plotting models b...")
    figs = plot_slices_temperature(
        models,
        AP.n_params,
        AP.temperature_bounds,
        AP.lattice_b_bounds,
        plot_bounds=[-10.0,320.0],
        property_name="Lattice b [angstrom]",
    )

    for fig in figs:
        pdf.savefig(fig)
    del figs

    # Plot parameter slices
    for param_name in AP.param_names:
        figs = plot_slices_params(
            models,
            param_name,
            AP.param_names,
            298,
            AP.temperature_bounds,
            AP.lattice_b_bounds,
            property_name="Lattice b [angstrom]",
        )
        for fig in figs:
            pdf.savefig(fig)
        del figs

    # Loop over test params
    for test_params in x_test[:,:AP.n_params]:
        train_points = []
        test_points = []
        # Locate rows where parameter set == test parameter set
        matches = np.unique(np.where((df_all[list(AP.param_names)] == test_params).all(axis=1))[0])
        # Loop over all matches -- these will be different temperatures
        for match in matches:
            # If the match (including T) is in the test set, then append to test points
            if np.where((df_all.values[match,:AP.n_params+1] == x_test[:,:AP.n_params+1]).all(axis=1))[0].shape[0] == 1:
                test_points.append([df_all["temperature"].iloc[match],df_all[property_name].iloc[match]])
            # Else append to train points
            else:
                train_points.append([df_all["temperature"].iloc[match],df_all[property_name].iloc[match]])

        pdf.savefig(
            plot_model_vs_test(
                models,
                test_params,
                np.asarray(train_points),
                np.asarray(test_points),
                AP.temperature_bounds,
                AP.lattice_b_bounds,
                plot_bounds=[-10.0, 320.0],
                property_name="Lattice b [angstrom]"
            )
        )

    ### Fit GP Model to lattice c
    param_names = list(AP.param_names) + ["temperature"]
    property_name = "sim_lattice_c"
    x_train, y_train, x_test, y_test = shuffle_and_split(
        df_all, param_names, property_name, shuffle_seed=gp_shuffle_seed, fraction_train=0.8
    )

    # Fit model
    print("Fitting models c")
    models = {}
    models["RBF"] = run_gpflow_scipy(
        x_train,
        y_train,
        gpflow.kernels.RBF(lengthscales=np.ones(AP.n_params + 1)),
    )
    models["Matern32"] = run_gpflow_scipy(
        x_train,
        y_train,
        gpflow.kernels.Matern32(lengthscales=np.ones(AP.n_params + 1)),
    )
    models["Matern52"] = run_gpflow_scipy(
        x_train,
        y_train,
        gpflow.kernels.Matern52(lengthscales=np.ones(AP.n_params + 1)),
    )

    # Plot model performance on train and test points
    pdf.savefig(plot_model_performance(models, x_train, y_train, AP.lattice_c_bounds))
    pdf.savefig(plot_model_performance(models, x_test, y_test, AP.lattice_c_bounds))

    # Plot temperature slices
    print("Plotting models c...")
    figs = plot_slices_temperature(
        models,
        AP.n_params,
        AP.temperature_bounds,
        AP.lattice_c_bounds,
        plot_bounds=[-10.0,320.0],
        property_name="Lattice c [angstrom]",
    )

    for fig in figs:
        pdf.savefig(fig)
    del figs

    # Plot parameter slices
    for param_name in AP.param_names:
        figs = plot_slices_params(
            models,
            param_name,
            AP.param_names,
            298,
            AP.temperature_bounds,
            AP.lattice_c_bounds,
            property_name="Lattice c [angstrom]",
        )
        for fig in figs:
            pdf.savefig(fig)
        del figs

    # Loop over test params
    for test_params in x_test[:,:AP.n_params]:
        train_points = []
        test_points = []
        # Locate rows where parameter set == test parameter set
        matches = np.unique(np.where((df_all[list(AP.param_names)] == test_params).all(axis=1))[0])
        # Loop over all matches -- these will be different temperatures
        for match in matches:
            # If the match (including T) is in the test set, then append to test points
            if np.where((df_all.values[match,:AP.n_params+1] == x_test[:,:AP.n_params+1]).all(axis=1))[0].shape[0] == 1:
                test_points.append([df_all["temperature"].iloc[match],df_all[property_name].iloc[match]])
            # Else append to train points
            else:
                train_points.append([df_all["temperature"].iloc[match],df_all[property_name].iloc[match]])

        pdf.savefig(
            plot_model_vs_test(
                models,
                test_params,
                np.asarray(train_points),
                np.asarray(test_points),
                AP.temperature_bounds,
                AP.lattice_c_bounds,
                plot_bounds=[-10.0, 320.0],
                property_name="Lattice c [angstrom]"
            )
        )


    pdf.close()


if __name__ == "__main__":
    main()
