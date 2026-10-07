from typing import List, Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def preparar_matriz_para_plot(matriz: np.ndarray, modo: str) -> np.ndarray:
    """Transforma la TDM al formato del modo de ploteo elegido."""
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
    tamano_ejes: int = 16,
    tamano_titulo: int = 14,
    guardar: bool = False,
    filename: Optional[str] = None,
    modo: str = "signed",
    color_grilla: str = "white",
    grosor_grilla: float = 1.5,
):
    """Heatmap en convención matemática estándar.

    * T[0,0] en la esquina superior-izquierda.
    * T[n-1,n-1] en la esquina inferior-derecha.
    * Diagonal principal: de arriba-izquierda a abajo-derecha.
    * Índice i (fila, donador) crece hacia ABAJO.
    * Índice j (columna, aceptor) crece hacia la DERECHA.
    * Ticks arriba (j) e izquierda (i).
    * Grilla entre celdas.
    """
    M_plot = preparar_matriz_para_plot(matriz, modo)

    # --- Rango de color ---
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

    n = M_plot.shape[0]
    if nombres is None:
        nombres = [str(i + 1) for i in range(n)]

    # --- Figura cuadrada ---
    fig, ax = plt.subplots(figsize=(6, 6))

    im = ax.imshow(
        M_plot,
        cmap=cmap,
        interpolation='nearest',
        vmin=vmin, vmax=vmax,
        origin='upper',              # ← T[0,0] arriba-izquierda
        aspect='equal',              # celdas cuadradas
    )

    # --- Ticks: arriba (j) e izquierda (i) ---
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(nombres, fontsize=tamano_ejes, fontweight='bold')
    ax.set_yticklabels(nombres, fontsize=tamano_ejes, fontweight='bold')

    ax.xaxis.tick_top()                          # ticks del eje j arriba
    ax.tick_params(axis='both', length=0)        # sin marquitas

    # --- Grilla entre celdas ---
    ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which='minor', color=color_grilla, linewidth=grosor_grilla)
    ax.tick_params(which='minor', length=0)

    # --- Colorbar ---
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    if titulo:
        ax.set_title(titulo, fontsize=tamano_titulo, pad=20)
    plt.tight_layout()

    if guardar and filename:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close(fig)