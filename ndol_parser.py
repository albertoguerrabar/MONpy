"""
Parser unificado para archivos de salida de NDOL 7.0.4 (Fortran).

Extrae:
  * Estructura molecular (coordenadas + número atómico + valencia)
  * Matriz completa de coeficientes MO (n_ao × n_mo)
  * Estados CI (singlete y triplete) con E, f_osc y transiciones c²_ia
  * Mapa AO → átomo

Maneja peculiaridades de Fortran: form feeds (\\x0c), \\r\\r\\n, \\x00, \\x1a.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Union

import numpy as np


# ----------------------------- Regex --------------------------------------

_RE_COORD_MARK = re.compile(r'INPUT ATOMIC CARTESIAN COORDINATES')
_RE_MO_MARK    = re.compile(r'SCF MOLECULAR ORBITALS')
_RE_COEF_TOKEN = re.compile(r'-?\d+\.\d+')
_RE_VALENCE    = re.compile(
    r'THIS SYSTEM HAS\s+(\d+)\s+VALENCE MOLECULAR ORBITALS,'
    r'\s+(\d+)\s+OF THEM OCCUPIED')

_RE_CI_SINGLET = re.compile(r'SINGLET\s*-\s*SINGLET TRANSITIONS')
_RE_CI_TRIPLET = re.compile(r'TRIPLET\s*-\s*TRIPLET TRANSITIONS')
_RE_CI_DETAIL  = re.compile(r'MAJOR ATOMIC ORBITAL CONTRIBUTIONS')
_RE_CI_HEADER  = re.compile(
    r'====>\s+SINGLY EXCITED CI STATE no\.\s+(\d+):\s+([\d.]+)\s+ev')
_RE_CI_TRANS   = re.compile(r"Single occ\. MO's\s+(\d+).*?=>\s+(\d+)")
_RE_CI_COEF    = re.compile(r'Sq\. CI coeff\.:\s+([\d.]+)')
_RE_CI_TABLE   = re.compile(
    r'^\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)')


# ----------------------------- Dataclasses --------------------------------

@dataclass
class Atomo:
    indice: int          # 0-based
    numero: int          # 1-based (columna ATOM)
    simbolo: str
    x: float
    y: float
    z: float
    nat: int
    valencia: int


@dataclass
class TransicionCI:
    i: int                       # MO ocupado (1-based)
    a: int                       # MO virtual (1-based)
    c2: float                    # Sq. CI coeff. = c_ia²

    @property
    def c_ia(self) -> float:
        """Amplitud CI (signo + por convención: NDOL no reporta el signo)."""
        return float(np.sqrt(max(self.c2, 0.0)))


@dataclass
class EstadoCI:
    indice: int
    energia_eV: float
    f_osc: float
    espin: str                    # 'singlet' | 'triplet'
    transiciones: List[TransicionCI] = field(default_factory=list)


@dataclass
class NDOLOutput:
    atomos: List[Atomo]
    C: np.ndarray                       # (n_ao, n_mo)
    energias_mo: np.ndarray             # (n_mo,)
    ao_a_atomo: np.ndarray              # (n_ao,) índices 0-based de átomo
    estados_singlete: List[EstadoCI] = field(default_factory=list)
    estados_triplete: List[EstadoCI] = field(default_factory=list)
    n_valencia_mo: Optional[int] = None
    n_ocupados: Optional[int] = None

    @property
    def n_ao(self) -> int:
        return self.C.shape[0]

    @property
    def n_mo(self) -> int:
        return self.C.shape[1]


# ----------------------------- Normalización Fortran ----------------------

def _read_ndol_lines(source: Union[str, bytes, List[str]]) -> List[str]:
    """Lee un archivo NDOL y normaliza peculiaridades de Fortran.

    * latin1 (ISO-8859-1) garantizado
    * Convierte form feeds (\\x0c) en saltos de línea
    * Elimina \\r, \\x00, \\x1a
    * Unifica saltos de línea
    """
    if isinstance(source, list):
        return source

    if isinstance(source, str):
        with open(source, 'rb') as f:
            raw = f.read()
    elif isinstance(source, bytes):
        raw = source
    else:
        raise TypeError(f"Tipo no soportado para source: {type(source)}")

    text = raw.decode('latin1', errors='replace')

    text = text.replace('\x0c', '\n')       # form feed → newline
    text = text.replace('\x1a', '')         # DOS EOF
    text = text.replace('\x00', '')         # NUL
    text = text.replace('\r\n', '\n').replace('\r', '\n')

    return text.split('\n')


# ----------------------------- Parsers internos ---------------------------

def _is_int_token(tok: str) -> bool:
    return bool(re.fullmatch(r'[+-]?\d+', tok))


def _parse_coordinates(lines: List[str]) -> List[Atomo]:
    i = 0
    while i < len(lines):
        if _RE_COORD_MARK.search(lines[i]):
            i += 2                      # saltar cabecera + columnas
            break
        i += 1
    else:
        raise ValueError("No se encontró bloque de coordenadas")

    atomos: List[Atomo] = []
    while i < len(lines):
        row = lines[i].strip()
        if not row:
            i += 1
            continue
        parts = row.split()
        if len(parts) < 6:
            break
        try:
            sym = parts[0]
            numero = int(parts[1])
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            nat = int(parts[5])
        except (ValueError, IndexError):
            break
        valencia = {'H': 1, 'C': 4, 'N': 4, 'O': 4}.get(sym, 0)
        atomos.append(Atomo(
            indice=len(atomos), numero=numero, simbolo=sym,
            x=x, y=y, z=z, nat=nat, valencia=valencia))
        i += 1
    return atomos


def _parse_mo_matrix(lines: List[str]) -> Tuple[np.ndarray, np.ndarray]:
    """Parser robusto de la sección SCF MOLECULAR ORBITALS.

    Estrategia:
      1. Detectar cabecera de bloque: línea con >=1 enteros y nada más.
      2. Saltar blancos, leer energías (mismo nº de floats que enteros).
      3. Saltar blancos, leer filas: cada una empieza por un entero
         (índice AO) seguido de n_cols floats.
      4. Continuar hasta agotar filas. Repetir con el siguiente bloque.
    """
    start = next(
        (i + 1 for i, ln in enumerate(lines) if _RE_MO_MARK.search(ln)),
        None,
    )
    if start is None:
        raise ValueError("'SCF MOLECULAR ORBITALS' no encontrado")

    # Localizar fin de la sección
    end = len(lines)
    for i in range(start, len(lines)):
        if 'CORE-VALENCE' in lines[i]:
            end = i
            break

    bloques_C: List[List[List[float]]] = []
    energias: List[float] = []

    i = start
    while i < end:
        raw = lines[i].strip()
        if not raw:
            i += 1
            continue

        tokens = raw.split()
        if not all(_is_int_token(t) for t in tokens):
            i += 1
            continue

        mo_idx = [int(t) for t in tokens]
        n_cols = len(mo_idx)

        # Buscar energías
        j = i + 1
        while j < end and not lines[j].strip():
            j += 1
        if j >= end:
            break
        try:
            E = [float(x) for x in lines[j].split()]
        except ValueError:
            i += 1
            continue
        if len(E) != n_cols:
            i += 1
            continue

        # Buscar filas de coeficientes
        j += 1
        while j < end and not lines[j].strip():
            j += 1

        bloque: List[List[float]] = []
        while j < end:
            row = lines[j].strip()
            if not row:
                j += 1
                continue

            parts = row.split(maxsplit=1)
            if not parts or not _is_int_token(parts[0]):
                break
            if len(parts) < 2:
                break

            coefs = [float(x) for x in _RE_COEF_TOKEN.findall(parts[1])]
            if len(coefs) != n_cols:
                break
            bloque.append(coefs)
            j += 1

        if bloque:
            bloques_C.append(bloque)
            energias.extend(E)
            i = j
        else:
            i += 1

    if not bloques_C:
        raise ValueError(
            "No se pudo parsear la matriz MO. "
            "Revisa el formato de la sección SCF MOLECULAR ORBITALS."
        )

    n_ao = len(bloques_C[0])
    for k, b in enumerate(bloques_C):
        if len(b) != n_ao:
            raise ValueError(
                f"Bloque {k}: {len(b)} filas != {n_ao} esperadas"
            )

    n_mo = len(energias)
    C = np.zeros((n_ao, n_mo), dtype=float)
    col = 0
    for bloque in bloques_C:
        n_cols = len(bloque[0])
        if col + n_cols > n_mo:
            raise ValueError(
                f"Exceso de columnas ({col + n_cols} > {n_mo})"
            )
        for r in range(n_ao):
            C[r, col:col + n_cols] = bloque[r]
        col += n_cols

    if col != n_mo:
        raise ValueError(f"Rellenadas {col} columnas vs {n_mo} energías")

    return C, np.array(energias)


def _build_ao_atom_map(atomos: List[Atomo]) -> np.ndarray:
    """Array (n_ao,) con el índice 0-based del átomo al que pertenece cada AO.

    Asume base mínima: H → 1 AO (s); resto → 4 AOs (s, px, py, pz).
    Si aparecen elementos con orbitales d, habría que extender.
    """
    mapping: List[int] = []
    for idx, at in enumerate(atomos):
        n = 1 if at.simbolo == 'H' else 4
        mapping.extend([idx] * n)
    return np.array(mapping, dtype=int)


def _parse_valence_info(lines: List[str]) -> Tuple[Optional[int], Optional[int]]:
    for line in lines:
        m = _RE_VALENCE.search(line)
        if m:
            return int(m.group(1)), int(m.group(2))
    return None, None


def _parse_ci_states(lines: List[str], spin: str) -> List[EstadoCI]:
    marker = _RE_CI_SINGLET if spin == 'singlet' else _RE_CI_TRIPLET

    start = None
    for i, line in enumerate(lines):
        if marker.search(line):
            start = i
            break
    if start is None:
        return []

    # Fin de la sección: siguiente *** CONFIGURATION INTERACTION ***
    section_end = len(lines)
    for j in range(start + 1, len(lines)):
        if '*** CONFIGURATION INTERACTION ***' in lines[j]:
            section_end = j
            break

    # 1) Tabla (E, f_osc)
    tabla = {}
    i = start + 1
    detail_start = None
    while i < section_end:
        if _RE_CI_DETAIL.search(lines[i]):
            detail_start = i
            break
        m = _RE_CI_TABLE.match(lines[i])
        if m:
            idx = int(m.group(1))
            E = float(m.group(2))
            f = float(m.group(6))
            tabla[idx] = (E, f)
        i += 1
    if detail_start is None:
        return []

    # 2) Bloque detallado
    estados: List[EstadoCI] = []
    j = detail_start + 1
    while j < section_end:
        m = _RE_CI_HEADER.search(lines[j])
        if not m:
            j += 1
            continue

        idx = int(m.group(1))
        E = float(m.group(2))
        f_osc = tabla.get(idx, (E, 0.0))[1]

        estado = EstadoCI(indice=idx, energia_eV=E, f_osc=f_osc, espin=spin)

        k = j + 1
        while k < section_end:
            if _RE_CI_HEADER.search(lines[k]):
                break
            mt = _RE_CI_TRANS.search(lines[k])
            if not mt:
                k += 1
                continue
            i_mo = int(mt.group(1))
            a_mo = int(mt.group(2))
            c2 = 0.0
            if k + 1 < section_end:
                mc = _RE_CI_COEF.search(lines[k + 1])
                if mc:
                    c2 = float(mc.group(1))
            estado.transiciones.append(TransicionCI(i=i_mo, a=a_mo, c2=c2))
            k += 1

        estados.append(estado)
        j = k

    return estados


# ----------------------------- API pública --------------------------------

def parse_ndol(source: Union[str, bytes, List[str]]) -> NDOLOutput:
    """Parsea un archivo NDOL (path, bytes, o lista de líneas)."""
    lines = _read_ndol_lines(source)

    atomos = _parse_coordinates(lines)
    C, E_mo = _parse_mo_matrix(lines)
    ao_a_atomo = _build_ao_atom_map(atomos)
    n_val, n_occ = _parse_valence_info(lines)
    estados_s = _parse_ci_states(lines, spin='singlet')
    estados_t = _parse_ci_states(lines, spin='triplet')

    if C.shape[0] != len(ao_a_atomo):
        raise ValueError(
            f"Inconsistencia: {C.shape[0]} AOs en C vs "
            f"{len(ao_a_atomo)} AOs esperados por los átomos "
            f"({len(atomos)} átomos)"
        )

    return NDOLOutput(
        atomos=atomos, C=C, energias_mo=E_mo, ao_a_atomo=ao_a_atomo,
        estados_singlete=estados_s, estados_triplete=estados_t,
        n_valencia_mo=n_val, n_ocupados=n_occ)