from typing import List, Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def preparar_matriz_para_plot(matriz: np.ndarray, modo: str) -> np.ndarray:
    """Transforma la TDM al formato del modo de ploteo elegido.

    Parameters
    ----------
    matriz : (n, n) ndarray
        TDM cruda (con signo)
    modo : str
        'signed'      → devuelve la matriz tal cual
        'magnitude'   → devuelve |T|
        'probability' → devuelve |T|² / Σ|T|²
    """
    if modo == "signed":
        return matriz.copy()
    if modo == "magnitude":
        return np.abs(matriz)
    if modo == "probability":
        M2 = matriz * matriz
        s = float(M2.sum())
        return M2 / s if s > 0 else M2
    raise ValueError(f"modo de ploteo desconocido: {modo}")


def graficar_matriz_densidad(
    matriz: np.ndarray,
    vmax_global: Optional[float] = None,
    nombres: Optional[List[str]] = None,
    titulo: str = "",
    cmap: str = "viridis",
    tamano_ejes: int = 18,
    tamano_titulo: int = 14,
    guardar: bool = False,
    filename: Optional[str] = None,
    modo: str = "signed",
):
    """Heatmap de la TDM reducida a fragmentos.

    Parameters
    ----------
    matriz : (n, n) ndarray
        TDM **cruda** (con signo). La transformación según `modo` se
        hace internamente.
    vmax_global : float o None
        Máximo de la escala de color. Para 'signed' el rango es
        [-vmax_global, +vmax_global]; para 'magnitude'/'probability'
        el rango es [0, vmax_global]. Si None, se calcula del propio
        `matriz`.
    modo : {'signed', 'magnitude', 'probability'}
        Cómo transformar y colorear la matriz.
    """
    M_plot = preparar_matriz_para_plot(matriz, modo)

    if modo == "signed":
        if vmax_global is None:
            vmax_global = float(np.abs(M_plot).max()) if M_plot.size else 1.0
        vmax_global = max(vmax_global, 1e-12)
        vmin, vmax = -vmax_global, vmax_global
    else:
        if vmax_global is None:
            vmax_global = float(M_plot.max()) if M_plot.size else 1.0
        vmax_global = max(vmax_global, 1e-12)
        vmin, vmax = 0.0, vmax_global

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(
        M_plot, cmap=cmap, interpolation='nearest',
        vmin=vmin, vmax=vmax,
    )

    n = M_plot.shape[0]
    if nombres is None:
        nombres = [str(i + 1) for i in range(n)]

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(nombres, fontsize=tamano_ejes, fontweight='bold')
    ax.set_yticklabels(nombres, fontsize=tamano_ejes, fontweight='bold')

    plt.colorbar(im, ax=ax)
    if titulo:
        ax.set_title(titulo, fontsize=tamano_titulo)
    plt.tight_layout()

    if guardar and filename:
        plt.savefig(filename, dpi=150)
    plt.close(fig)