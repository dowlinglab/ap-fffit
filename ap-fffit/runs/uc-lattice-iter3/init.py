import signac
import numpy as np
import pandas as pd
import unyt as u

from fffit.utils import values_scaled_to_real


def init_project():

    project = signac.init_project("ap-iter3")

    # Define temps
    temperatures = [
        298.0 * u.K,
        78.0 * u.K,
        10.0 * u.K,
    ]

    # Run at vapor pressure
    pressure = 1.0 * u.atm

    # Run for 200 ps (1 fs timestep)
    nsteps_eq = 100000
    nsteps_prod = 100000

    # Load samples from Latin hypercube
    lh_samples = pd.read_csv("/scratch365/bbefort/ap-fffit/ap-fffit/analysis/csv/uc-lattice-iter3-params.csv", header=0,usecols=[1,2,3,4,5,6,7,8]).values

    # Define bounds on sigma/epsilon
    # Sigma units = angstrom
    # Epsilon units = kcal/mol
    bounds_sigma = np.asarray(
        [
            [3.5, 4.5],  # Cl
            [0.5, 2.0],  # H
            [2.5, 3.8],  # N
            [2.5, 3.8],  # O
        ]
    )

    bounds_epsilon = np.asarray(
        [
            [0.1, 0.8],    # Cl
            [0.0, 0.02],   # H
            [0.01, 0.2],   # N
            [0.02, 0.3],   # O
        ]
    )

    bounds = np.vstack((bounds_sigma, bounds_epsilon))

    # Convert scaled latin hypercube samples to physical values
    scaled_params = values_scaled_to_real(lh_samples, bounds)

    for temperature in temperatures:
        for sample in scaled_params:

            (
                sigma_Cl,
                sigma_H,
                sigma_N,
                sigma_O,
                epsilon_Cl,
                epsilon_H,
                epsilon_N,
                epsilon_O,
            ) = sample

            state_point = {
                "T": float(temperature.to_value("K")),
                "P": float(pressure.to_value("atm")),
                "nsteps": {
                    "eq" : nsteps_eq,
                    "prod" : nsteps_prod,
                },
                "sigma_Cl": sigma_Cl,
                "sigma_H": sigma_H,
                "sigma_N": sigma_N,
                "sigma_O": sigma_O,
                "epsilon_Cl": epsilon_Cl,
                "epsilon_H": epsilon_H,
                "epsilon_N": epsilon_N,
                "epsilon_O": epsilon_O,
            }

            job = project.open_job(state_point)
            job.init()


if __name__ == "__main__":
    init_project()
