Installation
============

MONpy runs from its source directory. You need Python 3.10 or newer, NumPy and
Matplotlib. NDOL is required to produce a calculation output. Once you have the
``.ndo`` file, MONpy can analyze it independently.

Get the source
--------------

.. code-block:: bash

   git clone https://github.com/albertoguerrabar/MONpy.git
   cd MONpy

Create a Conda environment
--------------------------

If you already have a suitable environment, activate it and check that
``import numpy, matplotlib`` succeeds. To create a new one:

.. code-block:: bash

   conda create -n monpy -c conda-forge python=3.10 pip
   conda activate monpy
   python -m pip install -r requirements.txt

Run from a results directory
----------------------------

MONpy writes its results to your current directory. Create a directory for
each run to keep results from different calculations apart.

.. code-block:: bash

   monpy_dir="$(pwd)"
   mkdir -p results
   cd results
   python "$monpy_dir/main.py" /path/to/calculation.ndo

To run the bundled acrolein example from the repository root:

.. code-block:: bash

   python main.py --config examples/acrolein/monpy.toml

Omit the file argument to enter the path when MONpy starts. The program uses
Matplotlib's noninteractive ``Agg`` backend to save figures.
