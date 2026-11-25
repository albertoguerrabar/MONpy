"""
Parser unificado para archivos de salida de NDOL 7.0.4 (Fortran).

Extrae:
  * Estructura molecular (coordenadas + número atómico + valencia)
  * Matriz completa de coeficientes MO (n_ao × n_mo)
  * Estados CI (singlete y triplete) con E, f_osc y amplitudes CI cuando se imprimen
  * Mapa AO → átomo

Maneja peculiaridades de Fortran: form feeds (\\x0c), \\r\\r\\n, \\x00, \\x1a.
"""
from __future__ import annotations

import math
import re
from os import PathLike
from dataclasses import dataclass, field
from typing import Iterator, List, Optional, Tuple, Union

import numpy as np


# ----------------------------- Regex --------------------------------------

_RE_COORD_MARK = re.compile(r'INPUT ATOMIC CARTESIAN COORDINATES')
_RE_MO_MARK    = re.compile(r'SCF MOLECULAR ORBITALS')
_NUM = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
_RE_NUM = re.compile(_NUM)
_RE_COORD_HEADER = re.compile(r'^ATOM\s+X\s+Y\s+Z\s+NAT\b', re.I)
_RE_SYMBOL = re.compile(r'[A-Z][a-zA-Z]?')
_RE_VALENCE    = re.compile(
    r'THIS SYSTEM HAS\s+(\d+)\s+VALENCE MOLECULAR ORBITALS,'
    r'\s+(\d+)\s+OF THEM OCCUPIED')

_RE_CI_SINGLET = re.compile(r'SINGLET\s*-\s*SINGLET TRANSITIONS')
_RE_CI_TRIPLET = re.compile(r'TRIPLET\s*-\s*TRIPLET TRANSITIONS')
_RE_CI_DETAIL  = re.compile(r'MAJOR ATOMIC ORBITAL CONTRIBUTIONS')
_CI_HEADER = r'====>\s+SINGLY EXCITED CI STATE'
_RE_CI_HEADER_START = re.compile(_CI_HEADER, re.I)
_RE_CI_HEADER = re.compile(
    rf'{_CI_HEADER} no\.\s+(\d+):\s*({_NUM})\s+ev\b', re.I)
_CI_TRANS = r"Single occ\. MO's"
_RE_CI_TRANS_START = re.compile(_CI_TRANS)
_RE_CI_TRANS = re.compile(rf"{_CI_TRANS}\s+(\d+).*?=>\s+(\d+)(?=\s|$)")
_RE_CI_COEF = re.compile(rf'Sq\. CI coeff\.:\s*({_NUM})(?=\s|$)')
_RE_CI_FOSC = re.compile(rf'\bf\s*=\s*({_NUM})(?=\s|$)')
_RE_CI_SECTION = re.compile(r'\*\*\*\s+CONFIGURATION INTERACTION\s+\*\*\*')
_RE_CI_CAMPOS = re.compile(r'\s+|(?<=[\d.])(?=[+-])')
_RE_CI_VECTOR = re.compile(r'^\s*TRANSITION\s+CI COEFFICIENTS\s*$')
_RE_CI_VECTOR_ROW = re.compile(r'^\s*(\d+)\s+(\d+)\s*>\s*(\d+)\s+(.+)$')
_RE_CI_VECTOR_FIXED_ROW = re.compile(r'^(\*{3}|[ \d]{3})([ \d]{3})\s*>\s*(\d+)\s+(.+)$')
_RE_CI_N_CONFIG = re.compile(r'\b(\d+) CONFIGURATIONS ARE TAKEN INTO ACCOUNT\b')


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


@dataclass(slots=True)
class TransicionCI:
    i: int                       # MO ocupado (1-based)
    a: int                       # MO virtual (1-based)
    c2: float                    # Sq. CI coeff. = c_ia²
    amplitud: Optional[float] = None

    @property
    def c_ia(self) -> float:
        """Amplitud impresa o magnitud positiva si solo se conoce c_ia²."""
        if self.amplitud is not None:
            return self.amplitud
        return float(np.sqrt(max(self.c2, 0.0)))


@dataclass(slots=True)
class _TransicionesVectorCI:
    pares: tuple[tuple[int, int], ...]
    amplitudes: np.ndarray

    def __len__(self) -> int:
        return len(self.pares)

    def __iter__(self) -> Iterator[TransicionCI]:
        for (i, a), valor in zip(self.pares, self.amplitudes):
            amplitud = float(valor)
            yield TransicionCI(i, a, amplitud**2, amplitud)


@dataclass
class EstadoCI:
    indice: int
    energia_eV: float
    f_osc: float                 # NaN si NDOL no lo informa
    espin: str                    # 'singlet' | 'triplet'
    transiciones: List[TransicionCI] | _TransicionesVectorCI = field(default_factory=list)


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

def _read_ndol_lines(source: Union[str, PathLike[str], bytes, List[str]]) -> List[str]:
    """Lee un archivo NDOL y normaliza peculiaridades de Fortran.

    * latin1 (ISO-8859-1) garantizado
    * Convierte form feeds (\\x0c) en saltos de línea
    * Elimina \\r, \\x00, \\x1a
    * Unifica saltos de línea
    """
    if isinstance(source, list):
        text = '\n'.join(line.rstrip('\r\n') for line in source)
    elif isinstance(source, (str, PathLike)):
        with open(source, 'rb') as f:
            text = f.read().decode('latin1')
    elif isinstance(source, bytes):
        text = source.decode('latin1')
    else:
        raise TypeError(f"Tipo no soportado para source: {type(source)}")

    text = text.replace('\x0c', '\n')       # form feed → newline
    text = text.replace('\x1a', '')         # DOS EOF
    text = text.replace('\x00', '')         # NUL
    text = text.replace('\r\r\n', '\n').replace('\r\n', '\n').replace('\r', '\n')

    return text.split('\n')


# ----------------------------- Parsers internos ---------------------------

def _is_int_token(tok: str) -> bool:
    return bool(re.fullmatch(r'[+-]?\d+', tok))


def _parse_float(token: str, linea: int) -> float:
    if not _RE_NUM.fullmatch(token):
        raise ValueError(f"Línea {linea}: número inválido: {token!r}")
    valor = float(token.replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(valor):
        raise ValueError(f"Línea {linea}: número no finito: {token!r}")
    return valor


def _parse_numeric_row(texto: str, cantidad: int, linea: int) -> List[float]:
    valores = []
    fin = 0
    for match in _RE_NUM.finditer(texto):
        separador = texto[fin:match.start()]
        if separador.strip() or (
            valores and not separador and match[0][0] not in '+-'
        ):
            raise ValueError(f"Línea {linea}: campo numérico inválido: {texto!r}")
        valores.append(_parse_float(match[0], linea))
        fin = match.end()
    if texto[fin:].strip() or len(valores) != cantidad:
        raise ValueError(f"Línea {linea}: se esperaban {cantidad} campos numéricos")
    return valores


def _parse_coordinates(lines: List[str]) -> List[Atomo]:
    start = next((i + 1 for i, row in enumerate(lines)
                  if _RE_COORD_MARK.search(row)), None)
    if start is None:
        raise ValueError("No se encontró bloque de coordenadas")

    atomos: List[Atomo] = []
    numeros = set()
    for i in range(start, len(lines)):
        row = lines[i].strip()
        if not row or (not atomos and _RE_COORD_HEADER.match(row)):
            continue
        parts = row.split()
        if not _RE_SYMBOL.fullmatch(parts[0]):
            break
        if len(parts) < 6 or not all(_is_int_token(parts[j]) for j in (1, 5)):
            raise ValueError(f"Línea {i + 1}: fila de coordenadas inválida: {row!r}")
        sym = parts[0].capitalize()
        numero, nat = int(parts[1]), int(parts[5])
        if numero < 1 or numero in numeros:
            raise ValueError(f"Línea {i + 1}: número ATOM inválido o duplicado: {numero}")
        if nat < 1 or (sym == 'H') != (nat == 1):
            raise ValueError(f"Línea {i + 1}: NAT incompatible con {sym}: {nat}")
        x, y, z = (_parse_float(token, i + 1) for token in parts[2:5])
        valencia = 1 if nat == 1 else 4
        atomos.append(Atomo(len(atomos), numero, sym, x, y, z, nat, valencia))
        numeros.add(numero)
    if not atomos:
        raise ValueError("El bloque de coordenadas no contiene átomos")
    return atomos


def _parse_mo_matrix(
    lines: List[str], n_ao: int, n_mo_esperado: Optional[int]
) -> Tuple[np.ndarray, np.ndarray]:
    """Reconstruye C usando los números AO y MO impresos por NDOL."""
    start = next((i + 1 for i, row in enumerate(lines)
                  if _RE_MO_MARK.search(row)), None)
    if start is None:
        raise ValueError("'SCF MOLECULAR ORBITALS' no encontrado")

    end = next((i for i in range(start, len(lines))
                if 'CORE-VALENCE' in lines[i] or _RE_CI_SECTION.search(lines[i])),
               len(lines))
    filas = iter((i + 1, lines[i].strip()) for i in range(start, end)
                 if lines[i].strip())
    bloques = []
    energias = {}
    for linea, cabecera in filas:
        tokens = cabecera.split()
        if not all(_is_int_token(token) for token in tokens):
            if bloques:
                if _is_int_token(tokens[0]):
                    raise ValueError(f"Línea {linea}: fila AO fuera de un bloque MO")
                break
            continue
        mo = [int(token) for token in tokens]
        if len(set(mo)) != len(mo) or any(m < 1 or m in energias for m in mo):
            raise ValueError(f"Línea {linea}: índices MO inválidos o duplicados: {mo}")
        fila_energia = next(filas, None)
        if fila_energia is None:
            raise ValueError(f"Línea {linea}: bloque MO sin energías")
        numero_linea, texto = fila_energia
        energias.update(zip(mo, _parse_numeric_row(texto, len(mo), numero_linea)))

        bloque = np.empty((n_ao, len(mo)))
        vistos = set()
        for _ in range(n_ao):
            fila = next(filas, None)
            if fila is None:
                raise ValueError(f"Línea {linea}: bloque MO truncado, faltan filas AO")
            numero_linea, texto = fila
            parts = texto.split(maxsplit=1)
            if len(parts) != 2 or not _is_int_token(parts[0]):
                raise ValueError(f"Línea {numero_linea}: fila AO inválida: {texto!r}")
            ao = int(parts[0])
            if not 1 <= ao <= n_ao or ao in vistos:
                raise ValueError(f"Línea {numero_linea}: índice AO inválido o duplicado: {ao}")
            bloque[ao - 1] = _parse_numeric_row(parts[1], len(mo), numero_linea)
            vistos.add(ao)
        bloques.append((mo, bloque))

    if not bloques:
        raise ValueError("No se encontraron bloques de coeficientes MO")
    n_mo = len(energias) if n_mo_esperado is None else n_mo_esperado
    if len(energias) != n_mo or max(energias) != n_mo:
        raise ValueError(f"Se esperaban los índices MO consecutivos del 1 al {n_mo}")
    C = np.empty((n_ao, n_mo))
    for indices, bloque in bloques:
        C[:, np.asarray(indices) - 1] = bloque
    return C, np.array([energias[m] for m in range(1, n_mo + 1)])


def _build_ao_atom_map(atomos: List[Atomo]) -> np.ndarray:
    """Array (n_ao,) con el índice 0-based del átomo al que pertenece cada AO.

    Asume base mínima: H → 1 AO (s); resto → 4 AOs (s, px, py, pz).
    Si aparecen elementos con orbitales d, habría que extender.
    """
    mapping: List[int] = []
    for idx, at in enumerate(atomos):
        mapping.extend([idx] * at.valencia)
    return np.array(mapping, dtype=int)


def _parse_valence_info(lines: List[str]) -> Tuple[Optional[int], Optional[int]]:
    for line in lines:
        m = _RE_VALENCE.search(line)
        if m:
            return int(m.group(1)), int(m.group(2))
    return None, None


def _parse_ci_vectors(
    lines: List[str], start: int, end: int, spin: str,
    tabla: dict[int, EstadoCI], n_config: Optional[int],
) -> List[EstadoCI]:
    cabeceras = [j for j in range(start, end) if _RE_CI_VECTOR.match(lines[j])]
    estados: dict[int, EstadoCI] = {}
    configuraciones: Optional[dict[int, tuple[int, int]]] = None
    pares_ordenados: tuple[tuple[int, int], ...] = ()
    for principio, limite in zip(cabeceras, cabeceras[1:] + [end]):
        filas = iter((j + 1, lines[j].rstrip()) for j in range(principio + 1, limite)
                     if lines[j].strip())
        linea, cabecera = next(filas, (principio + 1, ''))
        tokens = cabecera.split()
        if not tokens or not all(token.isdigit() for token in tokens):
            raise ValueError(f"Línea {linea}: se esperaban números de estado CI")
        indices = [int(token) for token in tokens]
        if len(set(indices)) != len(indices) or any(i < 1 or i in estados for i in indices):
            raise ValueError(f"Línea {linea}: estados CI inválidos o duplicados: {indices}")

        linea, energias = next(filas, (linea + 1, ''))
        bloque = [
            EstadoCI(i, energia, tabla[i].f_osc if i in tabla else math.nan, spin)
            for i, energia in zip(indices, _parse_numeric_row(energias, len(indices), linea))
        ]
        pares = {}
        coeficientes = {}
        for linea, fila in filas:
            m = _RE_CI_VECTOR_ROW.match(fila) or _RE_CI_VECTOR_FIXED_ROW.match(fila)
            if m is None:
                if fila.split()[0].isdigit():
                    raise ValueError(f"Línea {linea}: fila de coeficientes CI inválida")
                break
            # La columna Fortran I3 desborda a *** a partir de la configuración 1000.
            # Solo se recupera por posición tras comprobar la secuencia anterior.
            if m[1] == "***":
                k = len(pares) + 1
                if k < 1000 or k - 1 not in pares:
                    raise ValueError(f"Línea {linea}: desbordamiento de índice CI inesperado")
            else:
                k = int(m[1])
            i_mo, a_mo = int(m[2]), int(m[3])
            if k < 1 or k in pares:
                raise ValueError(f"Línea {linea}: configuración CI inválida o duplicada: {k}")
            pares[k] = (i_mo, a_mo)
            valores = _parse_numeric_row(m[4], len(indices), linea)
            if any(abs(c) > 1 for c in valores):
                raise ValueError(f"Línea {linea}: amplitud CI fuera de [-1, 1]")
            coeficientes[k] = valores

        cantidad = len(pares) if n_config is None else n_config
        if not cantidad or pares.keys() != set(range(1, cantidad + 1)):
            raise ValueError(f"Línea {principio + 1}: faltan configuraciones en el bloque CI")
        if len(set(pares.values())) != cantidad:
            raise ValueError(f"Línea {principio + 1}: pares de MO duplicados en el bloque CI")
        if configuraciones is None:
            configuraciones = pares
            pares_ordenados = tuple(pares[k] for k in range(1, cantidad + 1))
        elif pares != configuraciones:
            raise ValueError(f"Línea {principio + 1}: configuraciones CI distintas entre bloques")

        amplitudes = np.array([coeficientes[k] for k in range(1, cantidad + 1)]).T.copy()
        for estado, vector in zip(bloque, amplitudes):
            estado.transiciones = _TransicionesVectorCI(pares_ordenados, vector)
        estados.update((estado.indice, estado) for estado in bloque)

    if tabla and estados.keys() != tabla.keys():
        raise ValueError(f"CI {spin}: los estados de CI COEFFICIENTS difieren de la tabla resumen")
    return [estados[i] for i in sorted(estados)]


def _parse_ci_states(
    lines: List[str], spin: str, n_config: Optional[int],
) -> List[EstadoCI]:
    marker = _RE_CI_SINGLET if spin == 'singlet' else _RE_CI_TRIPLET

    start = None
    for i, line in enumerate(lines):
        if marker.search(line):
            start = i
            break
    if start is None:
        return []

    section_end = next(
        (j for j in range(start + 1, len(lines))
         if (_RE_CI_SECTION.search(lines[j]) or _RE_CI_SINGLET.search(lines[j])
             or _RE_CI_TRIPLET.search(lines[j]))),
        len(lines),
    )

    tabla = {}
    columna_osc = None
    detail_start = section_end
    for i in range(start + 1, section_end):
        if (_RE_CI_DETAIL.search(lines[i]) or _RE_CI_HEADER_START.search(lines[i])
                or _RE_CI_VECTOR.match(lines[i])):
            detail_start = i
            break
        # Los campos Fortran pueden tocarse: «2.404-10.001» son dos valores.
        campos = _RE_CI_CAMPOS.split(lines[i].strip().upper())
        if campos[:2] == ['CI', 'E'] and 'OSC.' in campos:
            columna_osc = campos.index('OSC.')
        if columna_osc is None or not campos[0].isdigit():
            continue
        if len(campos) <= columna_osc:
            raise ValueError(f"Línea {i + 1}: fila de tabla CI incompleta")
        idx = int(campos[0])
        if idx in tabla:
            raise ValueError(f"Línea {i + 1}: estado CI {spin} duplicado: {idx}")
        tabla[idx] = EstadoCI(
            idx, _parse_float(campos[1], i + 1),
            _parse_float(campos[columna_osc], i + 1), spin,
        )

    if any(_RE_CI_VECTOR.match(lines[i]) for i in range(detail_start, section_end)):
        return _parse_ci_vectors(lines, detail_start, section_end, spin, tabla, n_config)

    estados = {}
    j = detail_start
    while j < section_end:
        m = _RE_CI_HEADER.search(lines[j])
        if not m:
            if _RE_CI_HEADER_START.search(lines[j]):
                raise ValueError(f"Línea {j + 1}: cabecera de estado CI inválida")
            j += 1
            continue
        idx = int(m[1])
        if idx in estados:
            raise ValueError(f"Línea {j + 1}: detalle CI {spin} duplicado: {idx}")
        f_osc = tabla[idx].f_osc if idx in tabla else math.nan
        estado = EstadoCI(idx, _parse_float(m[2], j + 1), f_osc, spin)
        estados[idx] = estado
        j += 1
        pares = set()
        while j < section_end and not _RE_CI_HEADER_START.search(lines[j]):
            fuerza = _RE_CI_FOSC.search(lines[j])
            if fuerza and idx not in tabla:
                estado.f_osc = _parse_float(fuerza[1], j + 1)
            transicion = _RE_CI_TRANS.search(lines[j])
            if not transicion:
                if _RE_CI_TRANS_START.search(lines[j]):
                    raise ValueError(f"Línea {j + 1}: transición CI inválida")
                j += 1
                continue
            i_mo, a_mo = int(transicion[1]), int(transicion[2])
            linea_transicion = j + 1
            j += 1
            while j < section_end and not lines[j].strip():
                j += 1
            coef = _RE_CI_COEF.search(lines[j]) if j < section_end else None
            if coef is None:
                raise ValueError(f"Línea {linea_transicion}: transición CI sin Sq. CI coeff.")
            c2 = _parse_float(coef[1], j + 1)
            if not 0 <= c2 <= 1:
                raise ValueError(f"Línea {j + 1}: coeficiente CI al cuadrado inválido: {c2}")
            if (i_mo, a_mo) in pares:
                raise ValueError(f"Línea {linea_transicion}: transición CI duplicada")
            pares.add((i_mo, a_mo))
            estado.transiciones.append(TransicionCI(i_mo, a_mo, c2))
            j += 1

    # Conserva la selección publicada por NDOL cuando hay algún detalle.
    return list((estados or tabla).values())


# ----------------------------- API pública --------------------------------

def parse_ndol(
    source: Union[str, PathLike[str], bytes, List[str]], *,
    espines: tuple[str, ...] = ('singlet', 'triplet'),
) -> NDOLOutput:
    """Parsea un archivo NDOL (path, bytes, o lista de líneas)."""
    lines = _read_ndol_lines(source)

    atomos = _parse_coordinates(lines)
    ao_a_atomo = _build_ao_atom_map(atomos)
    n_val, n_occ = _parse_valence_info(lines)
    C, E_mo = _parse_mo_matrix(lines, len(ao_a_atomo), n_val)
    if n_occ is not None and not 0 <= n_occ <= C.shape[1]:
        raise ValueError(f"Número inválido de MO ocupados: {n_occ}")

    n_config = next((int(m[1]) for line in lines if (m := _RE_CI_N_CONFIG.search(line))), None)
    estados_s = _parse_ci_states(lines, 'singlet', n_config) if 'singlet' in espines else []
    estados_t = _parse_ci_states(lines, 'triplet', n_config) if 'triplet' in espines else []
    pares_validados = set()
    for estado in estados_s + estados_t:
        transiciones = estado.transiciones
        if isinstance(transiciones, _TransicionesVectorCI):
            identificador = id(transiciones.pares)
            if identificador in pares_validados:
                continue
            pares_validados.add(identificador)
            pares = transiciones.pares
        else:
            pares = ((t.i, t.a) for t in transiciones)
        for i, a in pares:
            if (not 1 <= i < a <= C.shape[1]
                    or (n_occ is not None and not i <= n_occ < a)):
                raise ValueError(
                    f"CI {estado.espin} {estado.indice}: "
                    f"transición MO ocupado → virtual inválida: {i} → {a}"
                )

    return NDOLOutput(
        atomos=atomos, C=C, energias_mo=E_mo, ao_a_atomo=ao_a_atomo,
        estados_singlete=estados_s, estados_triplete=estados_t,
        n_valencia_mo=n_val, n_ocupados=n_occ)
