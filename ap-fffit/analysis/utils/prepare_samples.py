import numpy as np
import pandas as pd
from scipy.stats import linregress

from fffit.utils import values_scaled_to_real, values_real_to_scaled


def prepare_df_lattice(df_csv, molecule):
    """Prepare a pandas dataframe for fitting GP models to lattice constants

    Performs the following actions:
       - Removes any samples where simulation crashed
       - Converts all values from physical values to scaled values
       - Adds "expt_lattice_a"
       - Adds "expt_lattice_b"
       - Adds "expt_lattice_c"

    Parameters
    ----------
    df_csv : pd.DataFrame
        The dataframe as loaded from a CSV file with the signac results
    molecule : APConstants
        An instance of a molecule constants class
    n_molecules : int
        The number of molecules in the simulation

    Returns
    -------
    df_all : pd.DataFrame
        Dataframe with scaled parameters and simulation + expt. properties
        One row per simulation
    """
    if "a" not in df_csv.columns:
        raise ValueError("df_csv must contain column 'a'")
    if "b" not in df_csv.columns:
        raise ValueError("df_csv must contain column 'b'")
    if "c" not in df_csv.columns:
        raise ValueError("df_csv must contain column 'c'")
    for param in list(molecule.param_names):
        if param not in df_csv.columns:
            raise ValueError(
                f"df_csv must contain a column for parameter: '{param}'"
            )

    # Rename "lattice_xx" to "sim_lattice_xx"
    df_all = df_csv.rename(columns={
        "a": "sim_lattice_a",
        "b": "sim_lattice_b",
        "c": "sim_lattice_c",
        }
    )

    # Remove rows with missing data
    df_all = df_all.dropna()

    # Add the experimental results
    df_all["expt_lattice_a"] = df_all["temperature"].apply(
        lambda temp: molecule.expt_lattice_a[int(temp)]
    )
    df_all["expt_lattice_b"] = df_all["temperature"].apply(
        lambda temp: molecule.expt_lattice_b[int(temp)]
    )
    df_all["expt_lattice_c"] = df_all["temperature"].apply(
        lambda temp: molecule.expt_lattice_c[int(temp)]
    )

    # Scale ALL values
    scaled_param_values = values_real_to_scaled(
        df_all[list(molecule.param_names)], molecule.param_bounds
    )
    scaled_temperature = values_real_to_scaled(
        df_all["temperature"], molecule.temperature_bounds
    )
    scaled_sim_lattice_a = values_real_to_scaled(
        df_all["sim_lattice_a"], molecule.lattice_a_bounds
    )
    scaled_sim_lattice_b = values_real_to_scaled(
        df_all["sim_lattice_b"], molecule.lattice_b_bounds
    )
    scaled_sim_lattice_c = values_real_to_scaled(
        df_all["sim_lattice_c"], molecule.lattice_c_bounds
    )
    scaled_expt_lattice_a = values_real_to_scaled(
        df_all["expt_lattice_a"], molecule.lattice_a_bounds
    )
    scaled_expt_lattice_b = values_real_to_scaled(
        df_all["expt_lattice_b"], molecule.lattice_b_bounds
    )
    scaled_expt_lattice_c = values_real_to_scaled(
        df_all["expt_lattice_c"], molecule.lattice_c_bounds
    )


    # Construct the final df
    df_all[list(molecule.param_names)] = scaled_param_values
    df_all["temperature"] = scaled_temperature
    df_all["sim_lattice_a"] = scaled_sim_lattice_a
    df_all["sim_lattice_b"] = scaled_sim_lattice_b
    df_all["sim_lattice_c"] = scaled_sim_lattice_c
    df_all["expt_lattice_a"] = scaled_expt_lattice_a
    df_all["expt_lattice_b"] = scaled_expt_lattice_b
    df_all["expt_lattice_c"] = scaled_expt_lattice_c

    return df_all


def prepare_df_lattice_errors(df, molecule):
    """Create a dataframe with mean square error (mse) and mean absolute
    percent error (mape) for each unique parameter set.

    Parameters
    ----------
    df : pandas.Dataframe
        per simulation results
    molecule : APConstants
        molecule class with bounds and experimental data

    Returns
    -------
    df_new : pandas.Dataframe
        dataframe with one row per parameter set and
        the MSE and MAPE for each parameter set/property
    """
    new_data = []
    for group, values in df.groupby(list(molecule.param_names)):

        # Temperatures
        temps = values_scaled_to_real(
            values["temperature"], molecule.temperature_bounds
        )

        # Lattice a
        sim_lattice_a = values_scaled_to_real(
            values["sim_lattice_a"], molecule.lattice_a_bounds
        )
        expt_lattice_a = values_scaled_to_real(
            values["expt_lattice_a"], molecule.lattice_a_bounds
        )
        mse_lattice_a = np.mean((sim_lattice_a - expt_lattice_a) ** 2)
        mape_lattice_a = (
            np.mean(
                np.abs((sim_lattice_a - expt_lattice_a) / expt_lattice_a)
            )
            * 100.0
        )
        properties = {
            f"sim_lattice_a_{float(temp):.0f}K": float(lattice_a)
            for temp, lattice_a in zip(temps, sim_lattice_a)
        }

        # Lattice b
        sim_lattice_b = values_scaled_to_real(
            values["sim_lattice_b"], molecule.lattice_b_bounds
        )
        expt_lattice_b = values_scaled_to_real(
            values["expt_lattice_b"], molecule.lattice_b_bounds
        )
        mse_lattice_b = np.mean((sim_lattice_b - expt_lattice_b) ** 2)
        mape_lattice_b = (
            np.mean(
                np.abs((sim_lattice_b - expt_lattice_b) / expt_lattice_b)
            )
            * 100.0
        )
        properties.update(
            {
                f"sim_lattice_b_{float(temp):.0f}K": float(lattice_b)
                for temp, lattice_b in zip(temps, sim_lattice_b)
            }
        )

        # Lattice c
        sim_lattice_c = values_scaled_to_real(
            values["sim_lattice_c"], molecule.lattice_c_bounds
        )
        expt_lattice_c = values_scaled_to_real(
            values["expt_lattice_c"], molecule.lattice_c_bounds
        )
        mse_lattice_c = np.mean((sim_lattice_c - expt_lattice_c) ** 2)
        mape_lattice_c = (
            np.mean(
                np.abs((sim_lattice_c - expt_lattice_c) / expt_lattice_c)
            )
            * 100.0
        )
        properties.update(
            {
                f"sim_lattice_c_{float(temp):.0f}K": float(lattice_c)
                for temp, lattice_c in zip(temps, sim_lattice_c)
            }
        )

        new_quantities = {
            **properties,
            "mse_lattice_a": mse_lattice_a,
            "mse_lattice_b": mse_lattice_b,
            "mse_lattice_c": mse_lattice_c,
            "mape_lattice_a": mape_lattice_a,
            "mape_lattice_b": mape_lattice_b,
            "mape_lattice_c": mape_lattice_c,
        }

        new_data.append(list(group) + list(new_quantities.values()))

    columns = list(molecule.param_names) + list(new_quantities.keys())
    new_df = pd.DataFrame(new_data, columns=columns)

    # Remove parameter sets with missing data
    new_df = new_df.dropna()

    return new_df
