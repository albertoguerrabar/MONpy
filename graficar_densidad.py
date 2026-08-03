# graficar_densidad.py
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Optional

def plot_transition_density(matrix: np.ndarray,
                            vmax_global: Optional[float]=None,
                            names: Optional[List[str]]=None,
                            title: str="Transition Density Matrix",
                            cmap: str="viridis",
                            fontsize_ticks: int=18,
                            fontsize_title: int=14,
                            save: bool=False,
                            filename: Optional[str]=None):
    plt.figure(figsize=(8, 6))
    vmin = 0.0
    vmax = vmax_global if vmax_global is not None else np.max(matrix) if matrix.size > 0 else 1.0
    im = plt.imshow(matrix, cmap=cmap, interpolation='nearest', vmin=vmin, vmax=vmax)
    nrows, ncols = matrix.shape
    if names is None:
        names = [str(i+1) for i in range(ncols)]
    plt.xticks([])
    plt.yticks([])
    for idx, nm in enumerate(names):
        plt.text(idx, -0.6, nm, ha='center', va='bottom', fontsize=fontsize_ticks, fontweight='bold')
        plt.text(-0.6, idx, nm, ha='right', va='center', fontsize=fontsize_ticks, fontweight='bold')
    plt.colorbar(im, label="")
    plt.title(title, fontsize=fontsize_title)
    plt.tight_layout()
    if save and filename:
        plt.savefig(filename)
    plt.close()