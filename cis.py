import re
from typing import List, Dict, Any, Optional

# ------------------------------------------------------------------
# Regex precompiladas (una sola vez al importar el módulo)
# ------------------------------------------------------------------
_RE_CONFIG    = re.compile(r'\*\*\* CONFIGURATION INTERACTION \*\*\*')
_RE_SINGLET   = re.compile(r'SINGLET - SINGLET TRANSITIONS')
_RE_MAJOR     = re.compile(r'MAJOR ATOMIC ORBITAL CONTRIBUTIONS TO THE FIRST')
_RE_TRIPLET   = re.compile(r'TRIPLET - TRIPLET TRANSITIONS')
_RE_ESTADO    = re.compile(r'SINGLY EXCITED CI STATE no\.\s+(\d+):')
_RE_TRANS     = re.compile(r"Single occ\. MO's\s+(\d+).*=>\s+(\d+)")
_RE_COEF_CUAD = re.compile(r'Sq\. CI coeff\.\:\s+([0-9.]+)')
_RE_ATOMO     = re.compile(r'(\d+)[A-Z]')


def _parsear_float(s: str) -> float:
    s = s.strip()
    if not s:
        return 0.0
    try:
        return float(s.replace('D', 'E').replace('d', 'e'))
    except ValueError:
        return 0.0


def extraer_transiciones_con_coef(
    file_path: str,
    lines: Optional[List[str]] = None,
) -> List[Dict[str, Any]]:
    if lines is None:
        with open(file_path, 'r', encoding='latin1', errors='ignore') as f:
            lines = f.readlines()

    transiciones: Dict[int, List[Dict[str, Any]]] = {}
    estado_actual: Optional[int] = None
    dentro_config = dentro_singlet = dentro_major = False

    for i, line in enumerate(lines):
        if not dentro_config:
            if _RE_CONFIG.search(line):
                dentro_config = True
            continue

        if not dentro_singlet:
            if _RE_SINGLET.search(line):
                dentro_singlet = True
            continue

        if not dentro_major:
            if _RE_MAJOR.search(line):
                dentro_major = True
            continue

        if _RE_TRIPLET.search(line):
            break

        m_estado = _RE_ESTADO.search(line)
        if m_estado:
            estado_actual = int(m_estado.group(1))
            transiciones.setdefault(estado_actual, [])
            continue

        if estado_actual is None:
            continue

        m_trans = _RE_TRANS.search(line)
        if not m_trans:
            continue

        orbital_ocupado = int(m_trans.group(1))
        orbital_virtual = int(m_trans.group(2))

        # Coeficiente CI (c²) en la línea siguiente
        coef = 0.0
        if i + 1 < len(lines):
            m_coef = _RE_COEF_CUAD.search(lines[i + 1])
            if m_coef:
                coef = _parsear_float(m_coef.group(1))

        ocupados: List[tuple] = []
        virtuales: List[tuple] = []

        for j in range(2, 8):
            k = i + j
            if k >= len(lines):
                break
            row = lines[k]

            # --- OFFSETS IDÉNTICOS AL ORIGINAL (no tocar) ---
            m_at_occ = _RE_ATOMO.search(row[25:30])
            m_at_vir = _RE_ATOMO.search(row[53:58])

            atom_occ = int(m_at_occ.group(1)) if m_at_occ else None
            atom_vir = int(m_at_vir.group(1)) if m_at_vir else None

            coef_occ = _parsear_float(row[41:49])
            coef_vir = _parsear_float(row[69:77])

            ocupados.append((atom_occ, coef_occ))
            virtuales.append((atom_vir, coef_vir))

        transiciones[estado_actual].append({
            'orbital_ocupado': orbital_ocupado,
            'orbital_virtual': orbital_virtual,
            'coeficiente': coef,
            'ocupado': ocupados,
            'virtual': virtuales,
        })

    return [
        {"estado": eid, "transiciones": transiciones[eid]}
        for eid in sorted(transiciones)
    ]