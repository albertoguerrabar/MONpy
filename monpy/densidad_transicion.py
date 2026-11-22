"""Pesos electrón–hueco por bloques de la matriz de transición AO."""
from __future__ import annotations

import numpy as np

from .ndol_parser import EstadoCI, NDOLOutput


def construir_indices_ao_por_fragmento(
    ndol: NDOLOutput, fragmentos: dict[str, list[int]],
) -> list[np.ndarray]:
    """Convierte números ATOM de una partición molecular en índices de filas AO."""
    numeros = {atomo.numero for atomo in ndol.atomos}
    por_atomo: dict[int, int] = {}
    for indice, (nombre, atomos) in enumerate(fragmentos.items()):
        if not nombre.strip() or not atomos:
            raise ValueError("Cada fragmento debe tener nombre y al menos un átomo.")
        for numero in atomos:
            if numero not in numeros:
                raise ValueError(f"Átomo inexistente en {nombre}: {numero}")
            if numero in por_atomo:
                raise ValueError(f"Átomo repetido en la partición: {numero}")
            por_atomo[numero] = indice
    faltantes = numeros.difference(por_atomo)
    if faltantes:
        raise ValueError("Átomos sin fragmento: " + ", ".join(map(str, sorted(faltantes))))

    indices: list[list[int]] = [[] for _ in fragmentos]
    for ao, atomo in enumerate(ndol.ao_a_atomo):
        indices[por_atomo[ndol.atomos[atomo].numero]].append(ao)
    return [np.asarray(grupo, dtype=int) for grupo in indices]


def verificar_base_ortogonal(ndol: NDOLOutput) -> float:
    """Comprueba CᵀC ≈ I para la base efectiva impresa por NDOL."""
    if ndol.n_mo != ndol.n_ao:
        raise ValueError("Se requiere la matriz SCF cuadrada y completa.")
    error = float(np.max(np.abs(ndol.C.T @ ndol.C - np.eye(ndol.n_mo))))
    # Margen para el redondeo de los coeficientes SCF impresos.
    if not np.isfinite(error) or error > 1e-4:
        raise ValueError(f"CᵀC difiere de I ({error:.3g}). Se requiere información de solapamiento.")
    return error


def construir_tdm_ao(ndol: NDOLOutput, estado: EstadoCI) -> np.ndarray:
    """Construye D[AO del hueco, AO del electrón] a partir de amplitudes con signo."""
    if not estado.transiciones:
        raise ValueError("NDOL did not report CI transitions.")
    if any(t.amplitud is None for t in estado.transiciones):
        raise ValueError("Signed amplitudes missing. Enable CI COEFFICIENTS in NDOL.")

    i = np.array([t.i - 1 for t in estado.transiciones])
    a = np.array([t.a - 1 for t in estado.transiciones])
    if np.any((i < 0) | (i >= a) | (a >= ndol.n_mo)):
        raise ValueError("La expansión CI contiene índices MO inválidos.")

    # Comprime la expansión en los orbitales ocupados y virtuales que participan.
    ocupados, filas = np.unique(i, return_inverse=True)
    virtuales, columnas = np.unique(a, return_inverse=True)
    ci = np.zeros((ocupados.size, virtuales.size))
    np.add.at(ci, (filas, columnas), [t.c_ia for t in estado.transiciones])
    return ndol.C[:, ocupados] @ ci @ ndol.C[:, virtuales].T


def pesos_por_fragmento(D: np.ndarray, indices_ao: list[np.ndarray]) -> np.ndarray:
    """Suma los cuadrados de los elementos AO de cada bloque."""
    pesos_ao = np.square(D)
    return np.array([
        [pesos_ao[np.ix_(hueco, electron)].sum() for electron in indices_ao]
        for hueco in indices_ao
    ])


def analizar_pesos(omega: np.ndarray, nombres: list[str]) -> dict:
    total = float(omega.sum())
    if not np.isfinite(total):
        raise ValueError("El peso AO total no es finito.")
    if total == 0:
        return {"total": total, "fracciones": None}

    P = omega / total
    off = P.copy()
    np.fill_diagonal(off, 0)
    i, j = np.unravel_index(off.argmax(), off.shape)
    return {
        "total": total,
        "fracciones": P,
        "peso_LE": float(np.trace(P)),
        "Lambda_CT": float(off.sum()),
        "par_dominante": (nombres[i], nombres[j]) if off[i, j] > 0 else None,
        "peso_dominante": float(off[i, j]),
        "balance": P.sum(axis=0) - P.sum(axis=1),
    }
