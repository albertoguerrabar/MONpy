# propiedades.py
import numpy as np
from typing import Dict, Any

def exciton_size_aprox(M_norm: np.ndarray, positions: np.ndarray) -> float:
    """
    Calculate exciton size parameter d_exc (RMSdeh) - root-mean-square
    electron-hole distance using point-charge approximation.
    
    This is the physical electron-hole separation distance.
    
    Parameters:
    - M_norm: Normalized transition density matrix [n_frag, n_frag]
    - positions: Fragment centroid coordinates [n_frag, 3]
    
    Returns:
    - d_exc: RMS electron-hole distance in Angstroms
    """
    n = M_norm.shape[0]
    weighted_sum = 0.0
    
    for i in range(n):
        for j in range(n):
            # Distance between fragment i (hole) and fragment j (electron)
            r_vec = positions[i] - positions[j]
            distance_sq = np.dot(r_vec, r_vec)  # |r_i - r_j|²
            
            # Weight by normalized transition density probability
            weighted_sum += M_norm[i, j] * distance_sq
    
    # d_exc = sqrt(⟨|r_e - r_h|²⟩)
    d_exc = np.sqrt(weighted_sum)
    
    return float(d_exc)

def cti_absolute(M: np.ndarray) -> float:
    total = np.sum(np.abs(M))
    off = total - np.sum(np.abs(np.diag(M)))
    return 0.0 if total == 0 else float(off / total)

def cti_frobenius(M: np.ndarray) -> float:
    total = np.sum(M ** 2)
    off = total - np.sum(np.diag(M) ** 2)
    return 0.0 if total == 0 else float(off / total)

def fragment_contributions(M: np.ndarray):
    sq = M ** 2
    hole = np.sum(sq, axis=1)  # sum rows
    ele = np.sum(sq, axis=0)   # sum cols
    total = np.sum(sq)
    if total > 0:
        hole = hole / total
        ele = ele / total
    return hole, ele

def ifct_multiwfn(hole: np.ndarray, ele: np.ndarray) -> Dict[str, Any]:
    n = len(hole)
    if not np.isclose(np.sum(hole), 1.0, atol=1e-10):
        raise ValueError("Hole contributions do not sum to 1.")
    if not np.isclose(np.sum(ele), 1.0, atol=1e-10):
        raise ValueError("Electron contributions do not sum to 1.")

    Q = np.outer(hole, ele)
    p_net = np.zeros((n, n))
    for R in range(n):
        for S in range(n):
            if R != S:
                p_net[R, S] = Q[R, S] - Q[S, R]
    delta_p = np.zeros(n)
    for R in range(n):
        # population change = sum over S of p_{S->R} = sum_S p_net[S,R]
        delta_p[R] = np.sum([p_net[S, R] for S in range(n) if S != R])
    delta_verify = ele - hole
    if not np.allclose(delta_p, delta_verify, atol=1e-12):
        # warn but continue
        print("WARNING: IFCT verification mismatch (delta_p vs ele-hole).")

    intra = np.diag(Q)
    intrinsic_CT = 100.0 * (np.sum(Q) - np.sum(intra))
    intrinsic_LE = 100.0 * np.sum(intra)
    apparent_CT = None
    apparent_LE = None
    if n == 2:
        apparent_CT = 100.0 * abs(delta_p[0])
        apparent_LE = 100.0 - apparent_CT

    return {
        "Q_matrix": Q,
        "theta_hole": hole,
        "theta_ele": ele,
        "net_transfer": p_net,
        "population_change": delta_p,
        "population_change_verification": delta_verify,
        "intra_fragment": intra,
        "intrinsic_CT_percent": intrinsic_CT,
        "intrinsic_LE_percent": intrinsic_LE,
        "apparent_CT_percent": apparent_CT,
        "apparent_LE_percent": apparent_LE
    }

def interpret_exciton_size(dexc_value: float) -> str:
    """
    Interpret the exciton size based on d_exc value.
    """
    if dexc_value < 3.0:
        return "Localized (Frenkel) exciton"
    elif dexc_value < 10.0:
        return "Charge-transfer exciton"
    else:
        return "Separated (Wannier-Mott) exciton"

def calculate_all_properties_for_state(M_norm: np.ndarray, frag_centroids: np.ndarray) -> Dict[str, Any]:
    """
    Calculate all properties for an excited state using NORMALIZED matrix.
    
    Parameters:
    - M_norm: Normalized transition density matrix [n_frag, n_frag]
    - frag_centroids: Fragment centroid coordinates [n_frag, 3]
    
    Returns:
    - Dictionary with all computed properties
    """
    props = {}
    
    # Exciton size parameters
    props["exciton_size_aprox"] = exciton_size_aprox(M_norm, frag_centroids)
    props["exciton_size_interpretation"] = interpret_exciton_size(props["exciton_size_aprox"])
    
    # Charge transfer indices
    props["cti_absolute"] = cti_absolute(M_norm)
    props["cti_frobenius"] = cti_frobenius(M_norm)
    
    # Fragment contributions
    hole, ele = fragment_contributions(M_norm)
    props["fragment_hole"] = hole
    props["fragment_ele"] = ele
    props["ifct"] = ifct_multiwfn(hole, ele)
    
    return props