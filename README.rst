ap-fffit
==========

**ap-fffit** is a repository used for a case study of our
machine learning directed force field optimization workflow
outlined in `this preprint <https://arxiv.org/abs/2103.03208>`_.
Here we use the procedure to generate force fields for
Ammonium Perchlorate (AP).

Citation
~~~~~~~~
This work has been submitted for review. In the meantime, you
may cite the `preprint <https://arxiv.org/abs/2103.03208>`_ as:

.. code-block:: bash

    BJ Befort, RS DeFever, G Tow, AW Dowling, and EJ Maginn. Machine learning
    directed optimization of classical molecular modeling force fields. arXiv
    (2021), https://arxiv.org/abs/2103.03208


Available Data
~~~~~~~~~~~~~~

AP Parameter Sets
#################
The non-dominated and "top two" parameter sets for AP are
provided under ``ap-fffit/analysis/csv/``. The seventy parameter
sets that outperform our hand-tuned force field are found in
``uc-lattice-final-params.csv`` and the "top two" two non-dominated
parameter sets are found in ``ap-final-2.csv``. The parameter values
in the csv files are normalized between 0 and 1 based upon the
parameter bounds for each atom type (see manuscript, or
``ap-fffit/analysis/utils/ap.py`` for definitions).

Molecular Simulations
#####################
All molecular simulations were performed under ``ap-fffit/runs``.
Each iteration was managed with ``signac-flow``. Inside of each
directory in ``runs``, you will find all the necessary files to
run the simulations. Note that you may not get the exact same simulation
results due to differences in software versions, random seeds, etc.
Nonetheless, all of the results from our molecular simulations are saved
under ``ap-fffit/analysis/lattice-iterZZ-results.csv``, where
``ZZ`` is the iteration number.

Surrogate modeling
##################
All of the scripts for the surrogate modeling are provided in
``ap-fffit/analysis``, following the same naming structure as
the csv files.

Figures
#######
All scripts required to generate the primary figures in the
manuscript are reported under ``ap-fffit/final-figs`` and the
associated PDF files are located under ``ap-fffit/final-figs/pdfs``.

Using this package
~~~~~~~~~~~~~~~~~~

Installation
############

This package has a number of requirements that can be installed in
different ways. We recommend using a conda environment to manage
most of the installation and dependencies. However, some items will
need to be installed from source or pip.

We recommend starting with a fresh conda environment, then installing
the packages listed under ``requirements-pip.txt`` with pip, then
installing the packages listed under ``requirements-conda.txt`` with
conda, and finally installing a few items from source
``requirements-other.txt``. We recommend ``python3.7`` and
taking packages from ``conda-forge``.

Running the simulations will also require an installation of LAMMPS.
This can be installed separately (see installation instructions
`here <https://docs.lammps.org/Install.html>`_ )

**WARNING**: Cloning the ``ap-fffit`` repository will take some time
and ~1.5 GB of disk space since it contains the Latin hypercube
that have ~1e6 parameter sets each.

An example of the procedure is provided below:

.. code-block:: bash

    # First clone ap-fffit and install pip/conda available dependencies
    # with a new conda environment named ap-fffit
    git clone git@github.com:dowlinglab/ap-fffit.git
    cd ap-fffit/
    conda create --name ap-fffit python=3.7 -c conda-forge
    conda activate ap-fffit
    python3 -m pip install -r requirements-pip.txt
    conda install --file requirements-conda.txt -c conda-forge
    cd ../

    # Now clone and install  other dependencies
    git clone git@github.com:dowlinglab/fffit.git
    # Checkout the v0.1 release of fffit and install
    cd fffit/
    git checkout tags/v0.1
    pip install .
    cd ../
    # Checkout the v0.1 release of block average and install
    git clone git@github.com:rsdefever/block_average.git
    cd block_average/
    git checkout tags/v0.1
    pip install .
    cd ../
    # Checkout the v0.1 release of LAMMPS Thermo and install
    git clone git@github.com:rsdefever/lammps_thermo.git
    cd lammps_thermo/
    git checkout tags/v0.1
    pip install .
    cd ../

Compiling the analysis codes
#############################

There are two fortran analysis codes under ``ap-fffit/runs/codes`` that
must be compiled. You can use ``ifort`` or ``gfortran``. For example:

.. code-block:: bash

    ifort calc_htweaked_pc_uc.f90 -o calc_htweaked_pc_uc -no-wrap-margin
    ifort calc_hbond_hangle.f90 -o calc_hbond_hangle -no-wrap-margin

AP force field optimization
###########################

**NOTE**: We use signac and signac flow (`<https://signac.io/>`_)
to manage the setup and execution of the molecular simulations. These
instructions assume a working knowledge of that software.

The first iteration of the ammonium perchlorate simulations were
performed under the ``ap-fffit/runs/uc-lattice-iter1/``.
A Latin hypercube with 250 parameter sets exists under
``ap-fffit/runs/data/LHS_ap_iter1.csv``.
The signac workspace is created by ``ap-fffit/runs/uc-lattice-iter1/init.py``.

.. code-block:: bash

    cd ap-fffit/runs/uc-lattice-iter1/
    python init.py

The thermodynamic conditions for the simulations and the bounds for each parameter
(LJ sigma and epsilon for Cl, O, N, and H) are defined inside ``init.py``.

The simulation workflow is
defined in ``ap-fffit/runs/uc-lattice-iter1/project.py``. The flow operations
defined therein create the simulation input files, perform the simulations,
and run the analysis (calculating the lattice constants and unit cell mean distance).
In order to run these flow operations on a cluster with a job scheduler, it will be
necessary to edit the files under
``ap-fffit/runs/uc-lattice-iter1/templates/`` to be compatible with
your cluster. The signac documentation contains the necessary details.

Once the first iteration of simulations have completed (i.e., all the flow
operations are done), you can perform analysis. The necessary files are located
under ``ap-fffit/runs/analysis`` and ``ap-fffit/runs/analysis/uc-lattice-iter1``.
The first step is to extract the results from your signac project into a CSV file
so they can be stored and accessed more easily in the future. This step is
performed by ``extract.py``. The script requires the iteration number
as a command line argument.

**WARNING**: Running this script will overwrite your local copy of our simulation
results (stored as CSV files) with the results from your simulations.

To extract the results for iteration 1 run the following:

.. code-block:: bash

    cd ap-fffit/analysis/
    python extract.py 1


The CSV file with the results is saved under
``ap-fffit/analysis/csv/uc-lattice-iterXX-results.csv`` where ``XX``
is the iteration number.

The analysis is performed within a separate directory for each iteration.
For example, for the first iteration, it is performed under
``ap-fffit/analysis/uc-lattice-iter1``. The script ``id-new-samples.py``
loads the results from the CSV file, fits the SVM classifier and GP surrogate
models, loads the Latin hypercube with 1e6 prospective parameter sets,
and identifies the 250 new parameter sets to use for molecular simulations in
iteration 2. These parameter sets are saved to a CSV file:
``ap-fffit/analysis/csv/uc-lattice-iter2-params.csv``.

The second iteration of the liquid density simulations were
performed under the ``ap-fffit/runs/uc-lattice-iter2/``. The procedure
is the same as for iteration 1, but this time the force field parameters
are taken from: ``ap-fffit/analysis/csv/uc-lattice-iter2-params.csv``.
The procedure for analysis is likewise analogous to iteration 1, however,
note that in training the surrogate models,
``ap-fffit/runs/analysis/uc-lattice-iter2/id-new-samples.py`` now uses
the simulation results from both iterations 1 and 2.

Credits
~~~~~~~

This work was supported by the National Science Foundation
under grant NSF Award Number OAC-1835630 and NSF Award Number CBET-1917474
and the Air Force Office of Scientific Research under Contract
AFOSR FA9550-18-1-0321. Any opinions, findings, and conclusions
or recommendations expressed in this material are those of the
author(s) and do not necessarily reflect the views of the National
Science Foundation or Air Force Office of Scientific Research.
