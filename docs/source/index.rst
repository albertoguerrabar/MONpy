MONpy
=====

MONpy analyzes the electron–hole character of excited states from NDOL output
files (``.ndo``). It reconstructs the transition matrix in the atomic-orbital
(AO) basis from signed CI amplitudes and reports how much of each excitation is
local to a molecular fragment or involves charge transfer between fragments.

The analysis uses fragments you define by their printed NDOL atom numbers.
For a dimer, each molecule can be a fragment, even when its atom numbers are
interleaved with those of the other molecule.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   installation
   configuration
   outputs

.. toctree::
   :maxdepth: 2
   :caption: Reference

   TC
   citing
