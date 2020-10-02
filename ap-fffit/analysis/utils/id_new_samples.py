import numpy as np
import pandas as pd

from fffit.utils import values_real_to_scaled, values_scaled_to_real


def rank_samples(samples, gp_model, molecule, property_name):
    """Evalulate the GP model for a samples and return ranked results
    from lowest to highest MSE with experiment across the temperature range

    Parameters
    ----------
    samples : np.ndarray, shape=(n_samples, n_params)
        Samples to rank
    gp_model : gpflow.model
        GP model to predict the property_name of each sample
    molecule : R32Constants, R125Constants
        An instance of a molecule constants class
    property_name : string
        The name of the property of interest. Valid options are
        "lattice_a", "lattice_b", "lattice_c"

    Returns
    -------
    ranked_samples : pd.DataFrame
        Samples sorted by MSE
    """

    valid_property_names = [
        "sim_lattice_a",
        "sim_lattice_b",
        "sim_lattice_c",
    ]

    if property_name not in valid_property_names:
        raise ValueError(
            "Invalid property_name {}. Supported property_names are "
            "{}".format(property_name, valid_property_names)
        )

    temperature_bounds = molecule.temperature_bounds
    if property_name == "sim_lattice_a":
        expt_property = molecule.expt_lattice_a
        property_bounds = molecule.lattice_a_bounds
    elif property_name == "sim_lattice_b":
        expt_property = molecule.expt_lattice_b
        property_bounds = molecule.lattice_b_bounds
    elif property_name == "sim_lattice_c":
        expt_property = molecule.expt_lattice_c
        property_bounds = molecule.lattice_c_bounds

    # Apply GP model and calculate mean squared errors (MSE) between
    # GP model predictions and experimental data for all parameter samples
    mse = _calc_gp_mse(
        gp_model, samples, expt_property, property_bounds, temperature_bounds
    )
    # Make pandas dataframes, rank, and return
    samples_mse = np.hstack((samples, mse.reshape(-1, 1)))
    samples_mse = pd.DataFrame(
        samples_mse, columns=list(molecule.param_names) + ["mse"]
    )
    ranked_samples = samples_mse.sort_values("mse")

    return ranked_samples


def _calc_gp_mse(
    gp_model, samples, expt_property, property_bounds, temperature_bounds
):
    """Calculate the MSE between the GP model and experiment for samples"""

    all_errs = np.empty(shape=(samples.shape[0], len(expt_property.keys())))
    col_idx = 0
    for (temp, density) in expt_property.items():
        scaled_temp = values_real_to_scaled(temp, temperature_bounds)
        xx = np.hstack((samples, np.tile(scaled_temp, (samples.shape[0], 1))))
        means_scaled, vars_scaled = gp_model.predict_f(xx)
        means = values_scaled_to_real(means_scaled, property_bounds)
        err = means - density
        all_errs[:, col_idx] = err[:, 0]
        col_idx += 1

    return np.mean(all_errs ** 2, axis=1)
