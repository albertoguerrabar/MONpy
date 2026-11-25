from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, PowerNorm
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


def graficar_fracciones(
    fracciones: np.ndarray, nombres: list[str], filename: str | Path,
    *, color_gamma: float = 1.0,
) -> None:
    norm = (Normalize(vmin=0, vmax=1) if color_gamma == 1
            else PowerNorm(gamma=color_gamma, vmin=0, vmax=1))
    fig, ax = plt.subplots(figsize=(6, 6))
    try:
        imagen = ax.imshow(fracciones, cmap="viridis", interpolation="nearest",
                           norm=norm, origin="upper", aspect="equal")
        posiciones = range(len(nombres))
        ax.set_xticks(posiciones, nombres)
        ax.set_yticks(posiciones, nombres)
        ax.xaxis.tick_top()
        ax.tick_params(axis="both", which="major", labelsize=14, length=0, pad=7)

        limites = np.arange(len(nombres) - 1) + 0.5
        ax.set_xticks(limites, minor=True)
        ax.set_yticks(limites, minor=True)
        ax.tick_params(which="minor", length=0)
        ax.grid(which="minor", color="white", linewidth=0.5)

        eje_barra = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.14)
        fig.colorbar(imagen, cax=eje_barra)
        eje_barra.tick_params(labelsize=12)
        fig.savefig(filename, dpi=150, bbox_inches="tight")
    finally:
        plt.close(fig)
