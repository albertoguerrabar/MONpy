"""TDM reducida a fragmentos — versión corregida (fase + amplitudes CI).

T[F, G] = Σ_{i,a} c_ia · W_{F,i} · W_{G,a}

Esquemas de reducción disponibles:
  * 'signed'     : W_F,i = Σ_{μ∈F} C_{μi}      (con fase, riguroso)
  * 'population' : W_F,i = √(Σ_{μ∈F} C_{μi}²)  (Mulliken, énfasis CT)
  * 'abs_sum'    : W_F,i = Σ_{μ∈F} |C_{μi}|    (sin cancelación)
"""
from __future__ import annotations

from typing import Dict, List

import numpy as np

from ndol_parser import NDOLOutput, EstadoCI


# ---------------------------------------------------------------------------
# Construcción de índices
# ---------------------------------------------------------------------------

def construir_mapa_atomos_a_fragmento(
    fragmentos: Dict[str, List[int]]
) -> Dict[int, int]:
    """Número de átomo (1-based) → índice de fragmento (0-based)."""
    mapa: Dict[int, int] = {}
    for idx_f, (_, atomos) in enumerate(fragmentos.items()):
        for at in atomos:
            mapa[int(at)] = idx_f
    return mapa


def construir_indices_ao_por_fragmento(
    ndol: NDOLOutput,
    fragmentos: Dict[str, List[int]],
) -> List[np.ndarray]:
    """Precomputa los índices AO pertenecientes a cada fragmento."""
    n_frag = len(fragmentos)
    atom_to_frag = construir_mapa_atomos_a_fragmento(fragmentos)
    frag_ao: List[List[int]] = [[] for _ in range(n_frag)]

    for ao_idx, atom_idx in enumerate(ndol.ao_a_atomo):
        atom_num = int(atom_idx) + 1                # 1-based
        f = atom_to_frag.get(atom_num)
        if f is not None:
            frag_ao[f].append(ao_idx)

    return [np.array(x, dtype=int) for x in frag_ao]


# ---------------------------------------------------------------------------
# Reducción MO → fragmento según esquema
# ---------------------------------------------------------------------------

def _reducir_mo_a_fragmentos(
    C: np.ndarray,
    frag_ao_idx: List[np.ndarray],
    mo_idx: int,
    esquema: str,
) -> np.ndarray:
    """Devuelve W[F] para todos los fragmentos, según el esquema elegido.

    Parameters
    ----------
    C          : (n_ao, n_mo) matriz de coeficientes AO
    frag_ao_idx: lista de arrays con índices AO por fragmento
    mo_idx     : índice 0-based del MO a reducir
    esquema    : {'signed', 'population', 'abs_sum'}
    """
    n_frag = len(frag_ao_idx)
    W = np.zeros(n_frag, dtype=float)

    if esquema == 'signed':
        for k in range(n_frag):
            W[k] = C[frag_ao_idx[k], mo_idx].sum()
    elif esquema == 'population':
        for k in range(n_frag):
            W[k] = float(np.sqrt((C[frag_ao_idx[k], mo_idx] ** 2).sum()))
    elif esquema == 'abs_sum':
        for k in range(n_frag):
            W[k] = np.abs(C[frag_ao_idx[k], mo_idx]).sum()
    else:
        raise ValueError(f"esquema desconocido: {esquema}")

    return W


# ---------------------------------------------------------------------------
# TDM
# ---------------------------------------------------------------------------

def tdm_fragmentos(
    ndol: NDOLOutput,
    estado: EstadoCI,
    frag_ao_idx: List[np.ndarray],
    normalizar: str = 'none',
    esquema: str = 'signed',
) -> np.ndarray:
    """TDM (bloque occ→vir) reducida a fragmentos.

    Parameters
    ----------
    ndol        : NDOLOutput
    estado      : EstadoCI
    frag_ao_idx : índices AO por fragmento
    normalizar  : {'none', 'frobenius', 'suma_abs', 'max_abs', 'densidad'}
        Normalización numérica de la matriz resultante.
    esquema     : {'signed', 'population', 'abs_sum'}
        Cómo reducir cada MO al fragmento:
          * 'signed'     : W_F,i = Σ_{μ∈F} C_{μi}      (riguroso, con fase)
          * 'population' : W_F,i = √(Σ_{μ∈F} C_{μi}²)  (Mulliken)
          * 'abs_sum'    : W_F,i = Σ_{μ∈F} |C_{μi}|    (sin cancelación)

    Returns
    -------
    M : (n_frag, n_frag) ndarray
    """
    n_frag = len(frag_ao_idx)
    M = np.zeros((n_frag, n_frag), dtype=float)
    C = ndol.C

    for t in estado.transiciones:
        c_ia = t.c_ia
        if c_ia == 0.0:
            continue
        i_mo = t.i - 1
        a_mo = t.a - 1

        if not (0 <= i_mo < C.shape[1] and 0 <= a_mo < C.shape[1]):
            continue

        W_i = _reducir_mo_a_fragmentos(C, frag_ao_idx, i_mo, esquema)
        W_a = _reducir_mo_a_fragmentos(C, frag_ao_idx, a_mo, esquema)

        M += c_ia * np.outer(W_i, W_a)

    # --- Normalización numérica ---
    if normalizar == 'frobenius':
        nrm = float(np.sqrt((M * M).sum()))
        if nrm > 0:
            M = M / nrm
    elif normalizar == 'suma_abs':
        nrm = float(np.abs(M).sum())
        if nrm > 0:
            M = M / nrm
    elif normalizar == 'max_abs':
        nrm = float(np.abs(M).max()) if M.size else 0.0
        if nrm > 0:
            M = M / nrm
    elif normalizar == 'densidad':
        M2 = M * M
        nrm = float(M2.sum())
        if nrm > 0:
            M = M2 / nrm

    return M


# ---------------------------------------------------------------------------
# Análisis
# ---------------------------------------------------------------------------

def analizar_tdm(
    M: np.ndarray,
    nombres: List[str],
) -> Dict[str, object]:
    """Extrae métricas interpretativas de la TDM.

    Parameters
    ----------
    M       : matriz TDM (n_frag × n_frag)
    nombres : nombres de fragmentos

    Returns
    -------
    dict con:
        Lambda_CT     : índice de transferencia de carga ∈ [0,1]
        T_dominante   : valor del off-diagonal más grande en |·|
        par_dominante : (frag_sale, frag_llega)
        peso_LE       : Σ T_FF² / Σ T²  (peso local)
        dQ            : vector de ΔQ por fragmento (cols - filas)
        suma_filas    : S_F^out
        suma_cols     : S_G^in
    """
    M2 = M * M
    total = float(M2.sum())
    n = M.shape[0]

    if total == 0:
        return {
            'Lambda_CT': 0.0,
            'T_dominante': 0.0,
            'par_dominante': ('-', '-'),
            'peso_LE': 0.0,
            'dQ': np.zeros(n),
            'suma_filas': np.zeros(n),
            'suma_cols': np.zeros(n),
        }

    diag = float(np.trace(M2))
    lambda_ct = 1.0 - diag / total
    peso_le = diag / total

    # Par dominante off-diagonal
    off = M2 - np.diag(np.diag(M2))
    ij = np.unravel_index(int(np.argmax(off)), off.shape)
    T_dom = float(M[ij])
    par_dom = (nombres[ij[0]], nombres[ij[1]])

    suma_filas = M.sum(axis=1)     # cuánto sale de F
    suma_cols = M.sum(axis=0)      # cuánto llega a G
    dQ = suma_cols - suma_filas

    return {
        'Lambda_CT': lambda_ct,
        'T_dominante': T_dom,
        'par_dominante': par_dom,
        'peso_LE': peso_le,
        'dQ': dQ,
        'suma_filas': suma_filas,
        'suma_cols': suma_cols,
    }