ap-fffit
==========

**ap-fffit** is a repository used for a case study of our
machine learning directed force field optimization workflow
outlined in (this preprint). Here we use the procedure to generate
force fields for Ammonium Perchlorate (AP).

Citation
~~~~~~~~
This work has been submitted for review. In the meantime, you
may cite the `preprint <https://arxiv.org/abs/2103.03208>`_ as:

.. code-block:: bash

BJ Befort, RS DeFever, G Tow, AW Dowling, and EJ Maginn. Machine learning
directed optimization of classical molecular modeling force fields. arXiv
(2021), https://arxiv.org/abs/2103.03208


AP Parameter Sets
~~~~~~~~~~~~~~~~~
The non-dominated and "top two" parameter sets for AP are
provided under ``ap-fffit/analysis/csv/``. The seventy parameter
sets that outperform our hand-tuned force field are found in
``uc-lattice-final-params.csv`` and the "top two" two non-dominated 
parameter sets are found in ``ap-final-2.csv``. The parameter values
in the csv files are normalized between 0 and 1 based upon the
parameter bounds for each atom type (see manuscript, or
``ap-fffit/analysis/utils/ap.py`` for definitions).

Molecular Simulations
~~~~~~~~~~~~~~~~~~~~~
All molecular simulations were performed under ``ap-fffit/runs``.
Each iteration was managed with ``signac-flow``. Inside of each
directory in ``runs``, you will find all the necessary files to
run the simulations. Note that you may not get the exact same simulation
results due to differences in software versions, random seeds, etc.
Nonetheless, all of the results from our molecular simulations are saved
under ``ap-fffit/analysis/lattice-iterZZ-results.csv``, where
``ZZ`` is the iteration number.

Surrogate modeling
~~~~~~~~~~~~~~~~~~
All of the scripts for the surrogate modeling are provided in
``ap-fffit/analysis``, following the same naming structure as
the csv files.

Figures
~~~~~~~
All scripts required to generate the primary figures in the
manuscript are reported under ``ap-fffit/final-figs`` and the
associated PDF files are located under ``ap-fffit/final-figs/pdfs``.

Using this package
~~~~~~~~~~~~~~~~~~

This package has a variety of requirements that can be installed in
different ways. We recommend using a conda environment to manage
most of the installation and dependencies, and installing some items from
source or pip.

We recommend starting with a fresh conda environment, then installing
the packages listed under ``requirements-pip.txt`` with pip, then
installing the packages listed under ``requirements-conda.txt`` with
conda, and finally installing a few items from source
``requirements-other.txt``. We recommend ``python3.7`` and
taking packages from ``conda-forge``.

Credits
~~~~~~~

This work was supported by the National Science Foundation
under grant NSF Award Number OAC-1835630 and NSF Award Number CBET-1917474
and the Air Force Office of Scientific Research under Contract
AFOSR FA9550-18-1-0321. Any opinions, findings, and conclusions
or recommendations expressed in this material are those of the
author(s) and do not necessarily reflect the views of the National
Science Foundation or Air Force Office of Scientific Research.
