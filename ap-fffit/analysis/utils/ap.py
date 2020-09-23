import numpy as np
import unyt as u

class APConstants:
    """Experimental data and other constants for Ammonium Perchlorate"""
    def __init__(self):
        assert (
            self.expt_liq_density.keys()
            == self.expt_vap_density.keys()
            == self.expt_Pvap.keys()
            == self.expt_Hvap.keys()
        )

    @property
    def molecular_weight(self):
        """Molecular weight of the molecule in g/mol"""
        return 117.49

    @property
    def n_params(self):
        """Number of adjustable parameters"""
        return len(self.param_names)

    @property
    def param_names(self):
        """Adjustable parameter names"""

        param_names = (
            "sigma_Cl",
            "sigma_H",
            "sigma_N",
            "sigma_O",
            "epsilon_Cl",
            "epsilon_H",
            "epsilon_N",
            "epsilon_O",
        )

        return param_names

    @property
    def param_bounds(self):
        """Bounds on sigma and epsilon in units of nm and kJ/mol"""

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


    @property
    def expt_lattice_a(self):
        """Dictionary with experimental "a" lattice constant

        Temperature units K
        Lattice constant units Angstrom
        """

        expt_a = {
            10: 8.94,
            78: 9.02,
            298: 9.20,
        }

        return expt_a


    @property
    def expt_lattice_b(self):
        """Dictionary with experimental "b" lattice constant

        Temperature units K
        Lattice constant units Angstrom
        """

        expt_b = {
            10: 5.89,
            78: 5.85,
            298: 5.82,
        }

        return expt_b


    @property
    def expt_lattice_c(self):
        """Dictionary with experimental "c" lattice constant

        Temperature units K
        Lattice constant units Angstrom
        """

        expt_c = {
            10: 7.30,
            78: 7.39,
            298: 7.45,
        }

        return expt_c


    @property
    def temperature_bounds(self):
        """Bounds on temperature in units of K"""

        lower_bound = np.min(list(self.expt_lattice_a.keys()))
        upper_bound = np.max(list(self.expt_lattice_a.keys()))
        bounds = np.asarray([lower_bound, upper_bound], dtype=np.float32)
        return bounds

    @property
    def lattice_a_bounds(self):
        """Bounds on a lattice a constants in units of angstroms"""

        lower_bound = np.min(list(self.expt_lattice_a.values()))
        upper_bound = np.max(list(self.expt_lattice_a.values()))
        bounds = np.asarray([lower_bound, upper_bound], dtype=np.float32)
        return bounds

    @property
    def lattice_b_bounds(self):
        """Bounds on a lattice b constants in units of angstroms"""

        lower_bound = np.min(list(self.expt_lattice_b.values()))
        upper_bound = np.max(list(self.expt_lattice_b.values()))
        bounds = np.asarray([lower_bound, upper_bound], dtype=np.float32)
        return bounds

    @property
    def lattice_c_bounds(self):
        """Bounds on a lattice c constants in units of angstroms"""

        lower_bound = np.min(list(self.expt_lattice_c.values()))
        upper_bound = np.max(list(self.expt_lattice_c.values()))
        bounds = np.asarray([lower_bound, upper_bound], dtype=np.float32)
        return bounds

