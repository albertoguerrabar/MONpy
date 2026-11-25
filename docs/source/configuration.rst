Running an analysis
===================

Prepare the NDOL output
-----------------------

Use an ``.ndo`` file containing the Cartesian coordinates, the
``SCF MOLECULAR ORBITALS`` matrix and the ``CI COEFFICIENTS`` blocks for the
spin sections you want to analyze. The full CI blocks list the configurations
and their signed amplitudes. The abbreviated ``Sq. CI coeff.`` entries give
squared magnitudes. MONpy can read the abbreviated state information and
reports that quantitative CT is unavailable when the signed amplitudes are
missing.

Select states and fragments
---------------------------

#. Choose singlets, triplets or both.
#. Enter the displayed list positions separated by commas, or enter ``all``.
   The positions are distinct from the NDOL state numbers when both spins
   are listed together.
#. Enter the number and names of the fragments. For each fragment, use the
   numbers in the NDOL ``ATOM`` column. Commas and inclusive ranges are
   supported, such as ``1,3-6``. Every atom must belong to exactly one fragment.
#. Choose whether to export the signed AO transition matrices as ``.npy`` files.

Example with acrolein
---------------------

For an acrolein calculation with the atom order used below, an illustrative
partition groups atoms 1–4 and atoms 5–8. Use the printed atom numbers of your
own file when adapting the example.

.. code-block:: bash

   python /path/to/MONpy/main.py /path/to/acro.ndo

The following input selects the first singlet, labels the fragments A and B,
and uses the default of skipping AO matrix export.

.. code-block:: text

   Spin (1 singlet, 2 triplet, 3 both) [1]: 1
   Select states (e.g. 1,2 or all): 1
   Number of fragments: 2
   Fragment name 1: A
   ATOM numbers for 'A' (e.g. 1,3-6): 1-4
   Fragment name 2: B
   ATOM numbers for 'B' (e.g. 1,3-6): 5-8
   Export signed AO matrices (.npy)? [y/N]: n
   Color gamma (0.5 highlights small fractions) [1]: 1

Continue with :doc:`outputs` for the file names and matrix conventions.

Run from a TOML configuration
-----------------------------

Save these settings as ``monpy.toml`` alongside ``acro.ndo``.

.. code-block:: toml

   ndol = "acro.ndo"
   output_dir = "results/acro"
   export_ao = false
   color_gamma = 1.0

   [states]
   singlet = [1, 2]
   # triplet = "all"

   [fragments]
   A = "1-4"
   B = "5-8"

.. code-block:: bash

   python /path/to/MONpy/main.py --config monpy.toml

Paths in the TOML file are relative to its directory. State lists contain
NDOL state numbers. Set a spin to ``"all"`` to include every available state
of that spin, or omit it to skip that spin. Fragment order determines matrix
order. On Python 3.10, install ``tomli`` with ``python -m pip install tomli``.
Python 3.11 and newer include a TOML reader.

Use a positive ``color_gamma`` to control plot contrast. The default of 1
is linear, and 0.5 highlights small fractions. The numerical outputs use the
original fractions. Use the same gamma when comparing plots.
