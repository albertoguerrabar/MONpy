import sys

import numpy as np

from ndol_parser import parse_ndol
from fragmentos import definir_fragmentos
from densidad_transicion import (
    construir_indices_ao_por_fragmento,
    tdm_fragmentos,
    analizar_tdm,
)
from graficar_densidad import (
    graficar_matriz_densidad,
    preparar_matriz_para_plot,
)

LOGO = r"""
           __  __  ____  _   _
          |  \/  |/ __ \| \ | |
          | \  / | |  | |  \| |_ __  _   _
          | |\/| | |  | | . ` | '_ \| | | |
          | |  | | |__| | |\  | |_} | |_| |
          |_|  |_|\____/|_| \_| .__/ \__, |
                              | |     __/ |
                              |_|    |___/
"""

INFO_BLOCK = """
 *******************************************************
 *               Alberto Guerra-Barroso                *
 *                  Dept of Chemistry                  *
 *           Faculty of Exact Sciences, UNAB           *
 *                                                     *
 *                        and                          *
 *                                                     *
 *              Ana Lilian Montero-Alejo               *
 *            Dept of Physics, FCNMM, UTEM             *
 *******************************************************
"""


# ----------------------------------------------------------------------
# Utilidades de presentación
# ----------------------------------------------------------------------

def _formatear_matriz(matriz: np.ndarray, nombres) -> str:
    """Formatea la TDM como texto alineado, con cabecera de fragmentos."""
    header = " " * 12 + "".join(f"{n:>14}" for n in nombres)
    sep = "-" * (12 + 14 * len(nombres))
    filas = [header, sep]
    for i, nf in enumerate(nombres):
        fila = f"{nf:12}" + "".join(
            f"{matriz[i, j]:14.6f}" for j in range(len(nombres))
        )
        filas.append(fila)
    return "\n".join(filas)


def _formatear_analisis(info: dict, nombres: list) -> str:
    """Formatea el análisis de la TDM (Λ_CT, peso LE, par dominante, ΔQ)."""
    lineas = []
    lineas.append(f"  Lambda_CT = {info['Lambda_CT']:.4f}  "
                  f"(0 = LE puro, 1 = CT puro)")
    lineas.append(f"  peso LE   = {info['peso_LE']:.4f}")
    f1, f2 = info['par_dominante']
    lineas.append(f"  Par dominante: {f1} -> {f2}   "
                  f"T = {info['T_dominante']:+.4f}")
    lineas.append("  Delta Q por fragmento (heurístico, suma 0):")
    for k, nf in enumerate(nombres):
        lineas.append(f"      {nf:10s}  dQ = {info['dQ'][k]:+8.4f}")
    return "\n".join(lineas)


def _pedir_opcion(mensaje: str, opciones: dict, default: str) -> str:
    """Pide una opción mostrando el menú; devuelve la clave elegida."""
    print(mensaje)
    for k, v in opciones.items():
        print(f"{k} - {v}")
    raw = input(f"Select [{default}]: ").strip() or default
    return raw if raw in opciones else default


# ----------------------------------------------------------------------
# Flujo principal
# ----------------------------------------------------------------------

def main():
    print(LOGO)
    print("         MONpy: MO-based NDOL python\n")
    print(INFO_BLOCK)

    # --- Archivo de entrada ---
    archivo = (sys.argv[1] if len(sys.argv) > 1
               else input("Introduce the inputfile.ndo: ").strip())

    print(f"\nParsing {archivo}...")
    try:
        ndol = parse_ndol(archivo)
    except Exception as e:
        print(f"Error parsing NDOL output: {e}")
        return

    print(f"Molecular structure with {len(ndol.atomos)} atoms")
    print(f"MO matrix: {ndol.C.shape[0]} AOs x {ndol.C.shape[1]} MOs")
    if ndol.n_valencia_mo is not None:
        print(f"Total valence orbitals: {ndol.n_valencia_mo}, "
              f"{ndol.n_ocupados} occupied")
    print(f"Singlet CI states: {len(ndol.estados_singlete)}")
    print(f"Triplet CI states: {len(ndol.estados_triplete)}")

    # --- Fragmentos ---
    fragmentos = definir_fragmentos(ndol.atomos)
    if not fragmentos:
        print("Fragment mode is required.")
        return

    nombres_frag = list(fragmentos.keys())
    frag_ao_idx = construir_indices_ao_por_fragmento(ndol, fragmentos)

    print("\nFragments:")
    for nombre, atomos in fragmentos.items():
        print(f"  - {nombre}: {len(atomos)} atoms")

    # --- Elección de espín ---
    op_spin = _pedir_opcion(
        "\nSpin state:",
        {"1": "Singlet", "2": "Triplet", "3": "Both"},
        default="1",
    )
    if op_spin == "1":
        estados_todos = [(e, 'singlet') for e in ndol.estados_singlete]
    elif op_spin == "2":
        estados_todos = [(e, 'triplet') for e in ndol.estados_triplete]
    else:
        estados_todos = ([(e, 'singlet') for e in ndol.estados_singlete]
                         + [(e, 'triplet') for e in ndol.estados_triplete])

    if not estados_todos:
        print("No excited states found for the selected spin.")
        return

    # --- Selección de estados ---
    print(f"\nAvailable states: 1..{len(estados_todos)}")
    for k, (e, spin) in enumerate(estados_todos, start=1):
        print(f"  {k:3d}  [{spin:7s}]  E = {e.energia_eV:7.3f} eV   "
              f"f = {e.f_osc:6.4f}")

    entrada = input(
        "\nSelect states (e.g. 1,2 or 'all'): "
    ).strip().lower()
    if entrada == "all":
        indices = list(range(len(estados_todos)))
    else:
        try:
            indices = [int(i) - 1 for i in entrada.split(",") if i.strip()]
            indices = [i for i in indices if 0 <= i < len(estados_todos)]
            if not indices:
                print("No valid states selected.")
                return
        except ValueError:
            print("Invalid input.")
            return

    # --- Normalización numérica de la TDM ---
    op_norm = _pedir_opcion(
        "\nTDM numerical normalization:",
        {
            "1": "Frobenius   T / ||T||_F      (comparable entre estados)",
            "2": "MaxAbs      T / max|T|       (máximo = ±1, preserva signo)",
            "3": "SumaAbs     T / Σ|T|         (fracciones L1)",
            "4": "None        sin normalizar (raw)",
        },
        default="1",
    )
    norm_map = {
        "1": "frobenius",
        "2": "max_abs",
        "3": "suma_abs",
        "4": "none",
    }
    normalizacion = norm_map[op_norm]

    # --- Modos de ploteo a generar ---
    op_plot = _pedir_opcion(
        "\nPlot modes to generate:",
        {
            "1": "Signed only            (T, colorbar [-V, +V])",
            "2": "Magnitude only         (|T|, colorbar [0, 1])",
            "3": "Probability only       (|T|^2 / sum, colorbar [0, 1])",
            "4": "Both signed + magnitude",
            "5": "All three modes",
        },
        default="4",
    )
    modos = {
        "1": ["signed"],
        "2": ["magnitude"],
        "3": ["probability"],
        "4": ["signed", "magnitude"],
        "5": ["signed", "magnitude", "probability"],
    }[op_plot]

    # --- Cálculo de TDMs ---
    print("\nComputing transition density matrices...")
    resultados = []   # [(estado, spin, matriz | None, analisis | None)]

    for idx in indices:
        estado, spin = estados_todos[idx]
        etiqueta = f"{spin} #{estado.indice}"
        print(f"  {etiqueta}  (E = {estado.energia_eV:.3f} eV)...")
        try:
            M = tdm_fragmentos(
                ndol, estado, frag_ao_idx, normalizar=normalizacion,
            )
            info = analizar_tdm(M, nombres_frag)
            resultados.append((estado, spin, M, info))
        except Exception as e:
            print(f"    ERROR: {e}")
            resultados.append((estado, spin, None, None))

    # --- Vmax global por modo de ploteo ---
    matrices_validas = [M for _, _, M, _ in resultados if M is not None]

    vmax_por_modo = {}
    for modo in modos:
        if not matrices_validas:
            vmax_por_modo[modo] = 1.0
            continue
        transformed = [preparar_matriz_para_plot(M, modo)
                       for M in matrices_validas]
        if modo == "signed":
            vmax = max(float(np.abs(x).max()) for x in transformed)
        else:
            vmax = max(float(x.max()) for x in transformed)
        vmax_por_modo[modo] = max(vmax, 1e-12)

    print("\nGlobal color scales:")
    for modo, v in vmax_por_modo.items():
        if modo == "signed":
            print(f"  {modo:12s}  +/- {v:.6f}")
        else:
            print(f"  {modo:12s}  [0, {v:.6f}]")

    # --- Escritura ---
    output_filename = "monpy_results.txt"
    buf = [LOGO, "\n\n",
           "           MONpy: MO-based NDOL python\n\n",
           INFO_BLOCK, "\n\n"]

    buf.append(f"Input file: {archivo}\n")
    buf.append(f"Atoms: {len(ndol.atomos)}\n")
    buf.append(f"MO matrix: {ndol.C.shape[0]} AOs x {ndol.C.shape[1]} MOs\n")
    if ndol.n_valencia_mo is not None:
        buf.append(f"Valence orbitals: {ndol.n_valencia_mo} "
                   f"({ndol.n_ocupados} occupied)\n")
    buf.append(f"Fragments: {', '.join(nombres_frag)}\n")
    buf.append(f"Numerical normalization: {normalizacion}\n")
    buf.append(f"Plot modes: {', '.join(modos)}\n")
    buf.append("-" * 60 + "\n")

    for estado, spin, M, info in resultados:
        etiqueta = f"{spin.capitalize()} State {estado.indice}"
        header = (f"\n{etiqueta}:  "
                  f"E = {estado.energia_eV:.4f} eV, "
                  f"f = {estado.f_osc:.4f}\n")

        buf.append(header)
        print(header.strip())

        if M is None or info is None:
            buf.append("  (not computed)\n")
            continue

        # --- Análisis automático ---
        analisis_txt = _formatear_analisis(info, nombres_frag)
        buf.append("\nAnalysis:\n")
        buf.append(analisis_txt + "\n")
        print(analisis_txt)

        # --- Matriz completa en texto ---
        buf.append("\nTransition Density Matrix (numerical normalization "
                   f"= {normalizacion}):\n\n")
        buf.append(_formatear_matriz(M, nombres_frag))
        buf.append("\n" + "-" * 60 + "\n")

        print(np.round(M, 3))

        # --- Heatmaps (uno por modo) ---
        for modo in modos:
            nombre_png = f"TDM_{spin}_S{estado.indice}_{modo}.png"
            graficar_matriz_densidad(
                M,
                vmax_global=vmax_por_modo[modo],
                titulo="",
                nombres=nombres_frag,
                guardar=True,
                filename=nombre_png,
                modo=modo,
            )

    with open(output_filename, "w", encoding="utf-8") as f:
        f.write("".join(buf))

    print("\nHURRAY!!!")
    print(f"Results saved to {output_filename}")
    print("Images saved as TDM_<spin>_S<index>_<mode>.png")


if __name__ == "__main__":
    main()