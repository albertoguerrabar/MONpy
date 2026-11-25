import argparse
import csv
from pathlib import Path
import sys

import numpy as np

from .configuracion import leer_configuracion, seleccionar_estados
from .densidad_transicion import (
    analizar_pesos,
    construir_indices_ao_por_fragmento,
    construir_tdm_ao,
    pesos_por_fragmento,
    verificar_base_ortogonal,
)
from .fragmentos import definir_fragmentos, leer_numeros_atom
from .graficar_densidad import graficar_fracciones
from .ndol_parser import parse_ndol

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
 *    Lab of Computational and Theoretical Chemistry   *
 *      Faculty of Chemistry,University of Havana      *
 *                                                     *
 *                        and                          *
 *                                                     *
 *              Ana Lilian Montero-Alejo               *
 *          Department of Physics, FCNMM, UTEM         *
 *******************************************************
"""


CITATION_BLOCK = "Citation metadata is available in CITATION.cff."


def _formatear_matriz(matriz: np.ndarray, nombres: list[str]) -> str:
    ancho = max(14, max(map(len, nombres)) + 2)
    filas = [" " * ancho + "".join(f"{nombre:>{ancho}}" for nombre in nombres)]
    for nombre, fila in zip(nombres, matriz):
        filas.append(f"{nombre:<{ancho}}" + "".join(f"{valor:{ancho}.8g}" for valor in fila))
    return "\n".join(filas)


def _seleccionar_estados(ndol):
    opciones = {"1": {"singlet"}, "2": {"triplet"}, "3": {"singlet", "triplet"}}
    while True:
        opcion = input("Spin (1 singlet, 2 triplet, 3 both) [1]: ").strip() or "1"
        if opcion in opciones:
            break
        print("Choose 1, 2 or 3.")
    estados = [
        e for e in ndol.estados_singlete + ndol.estados_triplete
        if e.espin in opciones[opcion]
    ]
    for i, e in enumerate(estados, 1):
        f = f"{e.f_osc:.4f}" if np.isfinite(e.f_osc) else "unavailable"
        print(f"{i:4d}  {e.espin} #{e.indice}  E = {e.energia_eV:.4f} eV  f = {f}")
    while estados:
        seleccion = input("Select states (e.g. 1,2 or all): ").strip().lower()
        if seleccion == "all":
            return estados
        campos = [campo.strip() for campo in seleccion.split(",")]
        if all(campo.isdecimal() for campo in campos):
            indices = [int(campo) for campo in campos]
            if len(set(indices)) == len(indices) and all(1 <= i <= len(estados) for i in indices):
                return [estados[i - 1] for i in indices]
        print("Enter distinct numbers from the list, or all.")
    return []


def main() -> int:
    parser = argparse.ArgumentParser(description="Analyze NDOL excited states by fragment.")
    parser.add_argument("ndol", nargs="?", help="NDOL output file for interactive mode")
    parser.add_argument("--config", type=Path, help="TOML configuration for a noninteractive run")
    args = parser.parse_args()
    if args.config is not None and args.ndol is not None:
        parser.error("provide an NDOL output file or --config, not both")

    print(LOGO)
    print("         MONpy: MO-based NDOL python\n")
    print(INFO_BLOCK)
    try:
        config = leer_configuracion(args.config) if args.config is not None else None
        archivo = config.ndol if config else args.ndol or input("NDOL input file: ").strip()
        ndol = parse_ndol(archivo, espines=tuple(config.estados)) if config else parse_ndol(archivo)
    except (OSError, ValueError) as exc:
        print(f"Input error: {exc}", file=sys.stderr)
        return 1

    print(f"{len(ndol.atomos)} atoms | {ndol.n_ao} AOs × {ndol.n_mo} MOs")
    directorio = config.output_dir if config else Path.cwd()
    color_gamma = config.color_gamma if config else 1.0
    informe = [LOGO, INFO_BLOCK, CITATION_BLOCK, "", f"Input file: {archivo}",
               f"Atoms: {len(ndol.atomos)} | AOs: {ndol.n_ao} | MOs: {ndol.n_mo}"]
    try:
        if config:
            estados = seleccionar_estados(ndol, config)
            fragmentos = {}
            for nombre, entrada in config.fragmentos.items():
                try:
                    fragmentos[nombre] = leer_numeros_atom(entrada, len(ndol.atomos))
                except ValueError as exc:
                    raise ValueError(f"Fragment {nombre}: {exc}") from exc
            indices_ao = construir_indices_ao_por_fragmento(ndol, fragmentos)
            exportar_ao = config.export_ao
            directorio.mkdir(parents=True, exist_ok=True)
            informe.append(f"Configuration: {args.config.resolve()}")
        else:
            estados = _seleccionar_estados(ndol)
            fragmentos = {}
    except (OSError, ValueError) as exc:
        print(f"Could not prepare analysis: {exc}", file=sys.stderr)
        return 1

    if not estados:
        print("No states found for the selected spin.", file=sys.stderr)
        return 1

    con_amplitudes = [
        e for e in estados if e.transiciones
        and all(t.amplitud is not None for t in e.transiciones)
    ]
    error_base = None
    if con_amplitudes:
        try:
            error_ortogonalidad = verificar_base_ortogonal(ndol)
        except ValueError as exc:
            error_base = str(exc)
        else:
            if not config:
                fragmentos = definir_fragmentos(ndol.atomos)
                indices_ao = construir_indices_ao_por_fragmento(ndol, fragmentos)
                exportar_ao = input("Export signed AO matrices (.npy)? [y/N]: ").strip().lower() == "y"
                while True:
                    try:
                        gamma = float(input("Color gamma (0.5 highlights small fractions) [1]: ").strip() or "1")
                    except EOFError:
                        break
                    except ValueError:
                        gamma = 0
                    if np.isfinite(gamma) and gamma > 0:
                        color_gamma = gamma
                        break
                    print("Enter a finite positive number.")
            informe.append(f"Effective NDOL AO basis: max |CᵀC − I| = {error_ortogonalidad:.3e}")
            if color_gamma != 1.0:
                informe.append(f"Plot color gamma: {color_gamma:g}")
    nombres = list(fragmentos)
    if nombres:
        informe.append("Rows: hole fragment. Columns: electron fragment.")
        for i, (nombre, atomos) in enumerate(fragmentos.items(), 1):
            informe.append(f"Fragment {i} ({nombre}), ATOM: " + ", ".join(map(str, atomos)))

    filas_csv = []
    calculados = 0
    fallos_salida = False
    for estado in estados:
        etiqueta = f"{estado.espin.capitalize()} State {estado.indice}"
        f = f"{estado.f_osc:.4f}" if np.isfinite(estado.f_osc) else "unavailable"
        informe.append(f"\n{etiqueta} | E = {estado.energia_eV:.4f} eV | f = {f}")
        print(etiqueta)
        peso_ci = sum(t.c2 for t in estado.transiciones)
        fila = [estado.espin, estado.indice, estado.energia_eV,
                estado.f_osc if np.isfinite(estado.f_osc) else None,
                len(estado.transiciones), peso_ci]
        informe.append(f"  CI configurations printed: {len(estado.transiciones)} | sum(c^2) = {peso_ci:.6f}")
        if error_base or not estado.transiciones or any(t.amplitud is None for t in estado.transiciones):
            if not estado.transiciones:
                motivo = "NDOL did not report CI transitions."
            elif any(t.amplitud is None for t in estado.transiciones):
                motivo = "Signed amplitudes missing. Enable CI COEFFICIENTS in NDOL."
            else:
                motivo = error_base
            print(f"  {motivo}")
            informe.append(f"  CT unavailable. {motivo}")
            filas_csv.append(fila + [motivo, None, None, None, *([None] * len(nombres)**2)])
            continue

        D = construir_tdm_ao(ndol, estado)
        omega = pesos_por_fragmento(D, indices_ao)
        info = analizar_pesos(omega, nombres)
        informe.extend(["  Raw fragment weights (Omega)", _formatear_matriz(omega, nombres),
                        f"  Total weight: {info['total']:.8g}"])
        P = info["fracciones"]
        if P is None:
            motivo = "CT fractions undefined because the total weight is zero."
            print(f"  {motivo}")
            informe.append(f"  {motivo}")
            filas_csv.append(fila + [motivo, info["total"], None, None, *([None] * len(nombres)**2)])
        else:
            calculados += 1
            filas_csv.append(fila + ["ok", info["total"], info["Lambda_CT"], info["peso_LE"],
                                     *P.ravel().tolist()])
            metricas = f"  Lambda_CT = {info['Lambda_CT']:.6f} | LE = {info['peso_LE']:.6f}"
            print(metricas)
            informe.extend([metricas, "  Normalized fragment weights (P)", _formatear_matriz(P, nombres)])
            par = info["par_dominante"]
            if par is not None:
                informe.append(f"  Dominant CT pair: {par[0]} → {par[1]} ({info['peso_dominante']:.6f})")
            informe.append("  Electron–hole balance (fraction)")
            for nombre, balance in zip(nombres, info["balance"]):
                informe.append(f"    {nombre}: {balance:+.6f}")

        identificador = f"{estado.espin}_{estado.indice}"
        if exportar_ao:
            try:
                np.save(directorio / f"TDM_AO_{identificador}.npy", D)
            except OSError as exc:
                mensaje = f"Could not save AO matrix for {etiqueta}: {exc}"
                print(mensaje, file=sys.stderr)
                informe.append(mensaje)
                fallos_salida = True
        if P is not None:
            try:
                graficar_fracciones(
                    P, nombres, directorio / f"CT_{identificador}.png", color_gamma=color_gamma,
                )
            except OSError as exc:
                mensaje = f"Could not save CT plot for {etiqueta}: {exc}"
                print(mensaje, file=sys.stderr)
                informe.append(mensaje)
                fallos_salida = True

    try:
        with (directorio / "monpy_results.txt").open("w", encoding="utf-8") as salida:
            salida.write("\n".join(informe) + "\n")
        if config:
            with (directorio / "estados.csv").open("w", encoding="utf-8", newline="") as salida:
                writer = csv.writer(salida)
                writer.writerow([
                    "spin", "state", "E_eV", "f_osc", "ci_configurations", "sum_ci_squared",
                    "status", "omega_total", "lambda_CT", "peso_LE",
                    *(f"P[{i},{j}]" for i in range(1, len(nombres) + 1)
                      for j in range(1, len(nombres) + 1)),
                ])
                writer.writerows(filas_csv)
    except OSError as exc:
        print(f"Could not write results: {exc}", file=sys.stderr)
        return 1
    print(f"CT calculated for {calculados}/{len(estados)} states. See {directorio / 'monpy_results.txt'}.")
    return int(fallos_salida or not calculados)


if __name__ == "__main__":
    sys.exit(main())
