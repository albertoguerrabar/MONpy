# densidad_transicion.py
import numpy as np
from typing import Dict, List, Any

def build_fragment_transition_density(structure: List[Dict[str, Any]],
                                      fragments: Dict[int, List[int]],
                                      cis_data: Dict[str, Any],
                                      state_index: int = 0):
    """
    Build transition density matrix aggregated by fragments for a given excited state index (0-based).
    - structure: list of atom dicts (must contain 'number' for atom index)
    - fragments: dict mapping fragment_id -> list of atom indices (1-based)
    - cis_data: output from read_cis_transitions()
    - state_index: index in cis_data['excited_states'] list (0-based)
    Returns: matrix shape (n_frag, n_frag)
    """
    frag_names = list(fragments.keys())
    n_frag = len(frag_names)
    M = np.zeros((n_frag, n_frag), dtype=float)

    excited_states = cis_data.get("excited_states", [])
    if state_index < 0 or state_index >= len(excited_states):
        raise IndexError("state_index out of range in CIS data")

    state = excited_states[state_index]
    trans_list = state.get("transitions", [])

    # map atom -> fragment index for fast lookup
    atom_to_frag = {}
    for idx, frag in enumerate(frag_names):
        for atom in fragments[frag]:
            atom_to_frag[atom] = idx

    # Sum contributions: for each excitation contribution, add abs(coef * occ_coef * virt_coef)
    total_contrib = 0.0
    for trans in trans_list:
        coef = trans.get("coef", 0.0)
        occupied = trans.get("occupied", [])
        virtual = trans.get("virtual", [])
        # occupied and virtual are lists of (atom_index or None, coef)
        for occ_atom, occ_c in occupied:
            for virt_atom, virt_c in virtual:
                if occ_atom is None or virt_atom is None:
                    continue
                frag_i = atom_to_frag.get(occ_atom, None)
                frag_j = atom_to_frag.get(virt_atom, None)
                if frag_i is None or frag_j is None:
                    # atom not assigned to any fragment (shouldn't happen)
                    continue
                contrib = coef * (occ_c**2) * (virt_c**2)
                M[frag_i, frag_j] += contrib
                total_contrib += contrib

    # Normalización
    M_norm = np.zeros((n_frag, n_frag), dtype=float)
    if total_contrib > 0:
        M_norm = M / total_contrib
    else:
        M_norm = M
    
    return M, M_norm