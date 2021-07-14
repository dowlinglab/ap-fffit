## ap-fffit/runs

The directory naming scheme is ``uc-lattice-iter[iternum]`` where
``iternum`` is the iteration number.

Each directory contains the files necessary to run the molecular
simulations associated with the given iteration. The ``init.py``
defines the ``signac`` workspace and the ``project.py`` defines
the flow operations (i.e., performs the system setup, creates the
molecular simulation inputs, performs the simulations, and performs
any analysis).

