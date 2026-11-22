MONpy (MO-based NDOL Python) analyzes electron–hole localization and charge transfer from NDOL excited-state calculations. It reconstructs the transition matrix in the atomic-orbital (AO) basis from signed CI amplitudes, then computes local excitation and charge-transfer weights for user-defined molecular fragments.

## Quick start

MONpy requires Python 3.10 or newer, NumPy, and Matplotlib. Run it directly from the source checkout.

```bash
git clone https://github.com/albertoguerrabar/MONpy.git
cd MONpy
conda create -n monpy -c conda-forge python=3.10 numpy matplotlib
conda activate monpy

monpy_dir="$PWD"
mkdir -p results
cd results
python "$monpy_dir/main.py" /path/to/calculation.ndo
```

Use an NDOL output containing Cartesian coordinates, `SCF MOLECULAR ORBITALS`, and the `CI COEFFICIENTS` blocks for the spin sections of interest. The abbreviated `Sq. CI coeff.` entries can be read for state information. Quantitative charge-transfer analysis requires the signed amplitudes.

Select the spin and states, then define each fragment by the numbers printed in the NDOL `ATOM` column. Comma-separated numbers and inclusive ranges (for example, `1,3-6`) are supported. Assign each atom to exactly one fragment.

## Results

| File | Contents |
| --- | --- |
| `monpy_results.txt` | State energies and oscillator strengths, fragment definitions, raw and normalized fragment weights, and local excitation and charge-transfer descriptors. |
| `CT_<spin>_<state>.png` | Heatmap of normalized fragment weights, with a fixed scale from 0 to 1. |
| `TDM_AO_<spin>_<state>.npy` | Optional export of the signed AO transition matrix. |

Rows represent the hole fragment and columns the electron fragment. The diagonal gives the local excitation contributions, and the off-diagonal elements give interfragment charge transfer. The current implementation uses an effective orthonormal AO basis and checks the printed SCF coefficients for consistency with that convention.

## Citation

When reporting results obtained with MONpy, cite this repository and record the Git commit used for the calculation. The MONpy article reference and DOI are pending.
