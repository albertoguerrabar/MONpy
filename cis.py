# cis.py
import re
from typing import List, Dict, Any

def read_cis_transitions(file_path: str) -> Dict[str, Any]:
    """
    Parse CIS transitions from an .ndo-like file.
    Returns a dict with key "excited_states" that maps to a list of states.
    Each state is {"state": int, "transitions": [ ... ]}.
    """
    with open(file_path, 'r', encoding='latin1', errors='ignore') as f:
        lines = f.readlines()

    inside_config = False
    inside_singlet = False
    inside_major_block = False
    current_state = None
    transitions = {}
    i = 0
    while i < len(lines):
        line = lines[i]
        if '*** CONFIGURATION INTERACTION ***' in line:
            inside_config = True
        if inside_config and 'SINGLET - SINGLET TRANSITIONS' in line:
            inside_singlet = True
        if inside_singlet and 'MAJOR ATOMIC ORBITAL CONTRIBUTIONS TO THE FIRST' in line:
            inside_major_block = True

        if inside_major_block:
            m_state = re.search(r'SINGLY EXCITED CI STATE no\.\s+(\d+):', line)
            if m_state:
                current_state = int(m_state.group(1))
                transitions[current_state] = []
                i += 1
                continue

            m_trans = re.search(r"Single occ\. MO's\s+(\d+).*=>\s+(\d+)", line)
            if m_trans:
                occ_orb = int(m_trans.group(1))
                virt_orb = int(m_trans.group(2))
                coef = None
                occupied = []
                virtual = []

                # attempt to find coefficient on next line
                if i + 1 < len(lines):
                    m_coef = re.search(r'Sq\. CI coeff\.\:\s*([+-]?\d*\.\d+|[+-]?\d+)', lines[i + 1])
                    if m_coef:
                        try:
                            coef = float(m_coef.group(1))
                        except ValueError:
                            coef = None

                # parse following block lines for atomic contributions (best-effort)
                for j in range(2, 8):
                    if i + j < len(lines):
                        l = lines[i + j]
                        # try to extract atom numbers and coefficients at plausible columns
                        m_occ_atom = re.search(r'(\d+)[A-Z]', l[25:30])  # heuristic
                        m_virt_atom = re.search(r'(\d+)[A-Z]', l[53:58])
                        try:
                            occ_atom = int(m_occ_atom.group(1)) if m_occ_atom else None
                        except ValueError:
                            occ_atom = None
                        try:
                            virt_atom = int(m_virt_atom.group(1)) if m_virt_atom else None
                        except ValueError:
                            virt_atom = None

                        # coefficients heuristics (may vary by output)
                        try:
                            occ_coef = float(l[41:49].strip())
                        except Exception:
                            occ_coef = 0.0
                        try:
                            virt_coef = float(l[69:77].strip())
                        except Exception:
                            virt_coef = 0.0

                        occupied.append((occ_atom, occ_coef))
                        virtual.append((virt_atom, virt_coef))

                if coef is not None and current_state is not None:
                    transitions[current_state].append({
                        "occupied_orb": occ_orb,
                        "virtual_orb": virt_orb,
                        "coef": coef,
                        "occupied": occupied,   # list of (atom_index or None, coef)
                        "virtual": virtual
                    })

        if 'TRIPLET - TRIPLET TRANSITIONS' in line:
            break
        i += 1

    excited_states = []
    for st in sorted(transitions.keys()):
        excited_states.append({"state": st, "transitions": transitions[st]})

    return {"excited_states": excited_states}