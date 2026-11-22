Computational method
====================

Transition matrix
-----------------

For a selected state, MONpy places the signed CI amplitudes in a matrix
:math:`A_{\mathrm{CI}}` indexed by occupied and virtual molecular orbitals.
It then transforms the transition matrix to the AO representation:

.. math::

   D = C_{\mathrm{occ}} A_{\mathrm{CI}} C_{\mathrm{vir}}^{\mathsf{T}}.

Rows of :math:`D` represent the hole and columns represent the electron.
The current implementation uses the effective orthonormal AO convention and
checks :math:`C^{\mathsf{T}}C \approx I` against the printed SCF coefficients.
An explicit AO overlap correction is outside the current implementation.
The AO-to-atom map assumes the NDOL minimal valence basis, with one AO for H
and four for each other atom, grouped in coordinate order.

Fragment weights
----------------

For hole fragment :math:`F` and electron fragment :math:`G`, MONpy sums the
squared AO matrix elements of the corresponding block:

.. math::

   \Omega_{FG} = \sum_{\mu\in F,\,\nu\in G} |D_{\mu\nu}|^2,
   \qquad
   P_{FG} = \frac{\Omega_{FG}}{\sum_{F',G'}\Omega_{F'G'}}.

Here :math:`\Omega` contains the raw weights and :math:`P` contains their
normalized fractions. The fractions are undefined if the total weight is zero.
This analysis follows the fragment charge-transfer framework of
:ref:`Plasser and Lischka <method-references>` under the effective orthonormal
basis convention described above.

State descriptors
-----------------

* **Local excitation weight (LE)** is the sum of the diagonal elements of
  :math:`P`.
* **Charge-transfer weight** :math:`\Lambda_{\mathrm{CT}}` is the sum of
  its off-diagonal elements.
* **Dominant CT pair** is the largest positive off-diagonal element, reported
  from the hole fragment to the electron fragment.
* **Electron–hole balance** is the column sum minus the row sum for each
  fragment. It describes the difference of normalized marginals.

These descriptors characterize the CI expansion printed by NDOL. The
configuration count and :math:`\sum_{ia} c_{ia}^{2}` reported by MONpy provide
context about the expansion used.
