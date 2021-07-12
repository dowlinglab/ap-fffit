## ap-fffit/analysis

The directory naming scheme is ``uc-lattice-iter[iternum]`` where ``iternum`` is
the iteration number.

Each directory contains the files necessary to fit the support vector machine
(SVM) and Gaussian processs (GP) surrogate models, and then use the surrogate
models to identify the parameter sets to use for the next iteration of molecular
simulations.

The ``utils`` directory contains contains the following:

* ``ap.py``: parameter bounds and experimental reference data for each molecule
* ``save_signac_results.py``: helper functions to extract the simulation results from the signac projects
*  ``prepare_samples.py``: helper functions for preparing a ``pandas`` DataFrame with simulation results
*  ``id_new_samples.py``: helper functions for ranking the samples and computing the MSE for each parameter set
* ``plot.py``: helper functions for creating plots

The ``uc-lattice-final-analysis`` directory contains a script that extracts the
top-performing parameter sets, the ``csv`` directory contains CSV files
storing parameter sets and simulation results for each iteration, and the
``final-figs`` directory contains scripts and PDFs for the AP-related
figures found in the manuscript.

