import numpy as np
import unyt as u


class APConstants:
    """Experimental data and other constants for Ammonium Perchlorate"""
    def __init__(self):
        assert (
            self.expt_lattice_a.keys()
            == self.expt_lattice_b.keys()
            == self.expt_lattice_c.keys()
        )
        assert(
            len(self.primitive_cell_elements)
            == self.primitive_cell_xyz.shape[0]
        )
        assert(
            len(self.unit_cell_elements)
            == self.unit_cell_xyz.shape[0]
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

        bounds = np.vstack((bounds_sigma, bounds_epsilon))
        return bounds


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
    def primitive_cell_elements(self):
        elements = np.array([ "Cl", "O", "O", "O", "N", "H", "H", "H" ])
        return elements

    @property
    def unit_cell_elements(self):
        elements = np.array(
                [ "Cl", "O", "O", "O", "N", "H", "H", "H",
                  "Cl", "O", "O", "O", "N", "H", "H", "H",
                  "O", "H",
                  "Cl", "O", "O", "O", "N", "H", "H", "H",
                  "Cl", "O", "O", "O", "N", "H", "H", "H",
                  "O", "H", "O", "H", "O", "H",
                ]
        )
        return elements

    @property
    def primitive_cell_xyz(self):
        xyz = np.array(
            [
              [3.79592400000000, 1.47250000000000, 1.40963000000000],
              [2.62836000000000, 1.47250000000000, 0.55188000000000],
              [4.99388400000000, 1.47250000000000, 0.61028000000000],
              [3.76463400000000, 0.29685600000000, 2.25278000000000],
              [2.84739000000000, 1.47250000000000, 4.90998000000000],
              [1.90770909169815, 1.47250000000000, 4.49310617653451],
              [2.73888458499779, 1.47250000000000, 5.96240125991545],
              [3.38685404561816, 2.29027755286932, 4.58876441060340],
            ]
        )
        return xyz

    @property
    def unit_cell_xyz(self):
        xyz = np.array(
            [
                [3.79592400000000, 1.47250000000000, 1.40963000000000],
                [2.62836000000000, 1.47250000000000, 0.55188000000000],
                [4.99388400000000, 1.47250000000000, 0.61028000000000],
                [3.76463400000000, 0.29685600000000, 2.25278000000000],
                [2.84739000000000, 1.47250000000000, 4.90998000000000],
                [1.90770909169815, 1.47250000000000, 4.49310617653451],
                [2.73888458499779, 1.47250000000000, 5.96240125991545],
                [3.38685404561816, 2.29027755286932, 4.58876441060340],
                [0.67407600000000, 4.41750000000000, 5.05963000000000],
                [1.84164000000000, 4.41750000000000, 4.20188000000000],
                [8.41611600000000, 4.41750000000000, 4.26028000000000],
                [0.70536600000000, 3.24185600000000, 5.90278000000000],
                [1.62261000000000, 4.41750000000000, 1.25998000000000],
                [2.56229090830185, 4.41750000000000, 0.84310617653451],
                [1.73111541500221, 4.41750000000000, 2.31240125991545],
                [1.08314595438184, 5.23527755286932, 0.93876441060340],
                [3.76463400000000, 2.64814400000000, 2.25278000000000],
                [3.38685404561816, 0.65472244713068, 4.58876441060340],
                [8.26592400000000, 1.47250000000000, 2.24037000000000],
                [7.09836000000000, 1.47250000000000, 3.09812000000000],
                [0.52388400000000, 1.47250000000000, 3.03972000000000],
                [8.23463400000000, 0.29685600000000, 1.39722000000000],
                [7.31739000000000, 1.47250000000000, 6.04002000000000],
                [6.37770909169815, 1.47250000000000, 6.45689382346549],
                [7.20888458499779, 1.47250000000000, 4.98759874008455],
                [7.85685404561816, 2.29027755286932, 6.36123558939660],
                [5.14407600000000, 4.41750000000000, 5.89037000000000],
                [6.31164000000000, 4.41750000000000, 6.74812000000000],
                [3.94611600000000, 4.41750000000000, 6.68972000000000],
                [5.17536600000000, 5.59314400000000, 5.04722000000000],
                [6.09261000000000, 4.41750000000000, 2.39002000000000],
                [7.03229090830185, 4.41750000000000, 2.80689382346549],
                [6.20111541500221, 4.41750000000000, 1.33759874008455],
                [5.55314595438184, 3.59972244713068, 2.71123558939660],
                [8.23463400000000, 2.64814400000000, 1.39722000000000],
                [7.85685404561816, 0.65472244713068, 6.36123558939660],
                [5.17536600000000, 3.24185600000000, 5.04722000000000],
                [5.55314595438184, 5.23527755286932, 2.71123558939660],
                [0.70536600000000, 5.59314400000000, 5.90278000000000],
                [1.08314595438184, 3.59972244713068, 0.93876441060340],
            ]
        )

        return xyz

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

AP = APConstants()
