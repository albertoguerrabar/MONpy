import re
from typing import List, Dict, Any, Optional, Tuple

# Basic valence table (can be extended as needed)
VALENCE_TABLE = {
    "H": 1, "He": 0,
    "C": 4, "N": 3, "O": 2, "F": 1,
    "P": 3, "S": 2,
    "Cl": 1, "Br": 1, "I": 1,
}

# Regex for atomic coordinate lines:
COORD_REGEX = re.compile(
    r"^\s*([A-Z][a-z]?)\s+(\d+)\s+([+-]?\d*\.\d+|[+-]?\d+)\s+([+-]?\d*\.\d+|[+-]?\d+)\s+([+-]?\d*\.\d+|[+-]?\d+)\s+(\d+)"
)


def read_molecular_structure(lines: List[str]) -> List[Dict[str, Any]]:
    """
    Reads the Cartesian coordinate block from the .ndo file.
    Returns the molecular structure as a list of atom dictionaries.
    """
    structure = []
    inside_coords = False
    start_idx = None

    # Locate the beginning of the coordinate block
    for i, line in enumerate(lines):
        if "INPUT ATOMIC CARTESIAN COORDINATES" in line.upper():
            inside_coords = True
            continue

        if inside_coords and "ATOM" in line and "X" in line and "Y" in line and "Z" in line:
            start_idx = i + 2
            break

    if start_idx is None:
        raise ValueError("Coordinate block not found in .ndo file.")

    # Read coordinate lines
    for line in lines[start_idx:]:
        stripped = line.strip()
        if stripped == "":
            break

        match = COORD_REGEX.match(line)
        if not match:
            break

        symbol, number, x, y, z, nat = match.groups()

        structure.append({
            "symbol": symbol,
            "number": int(number),
            "x": float(x),
            "y": float(y),
            "z": float(z),
            "nat": int(nat),
            "valence": VALENCE_TABLE.get(symbol, 0)
        })

    return structure


# Multiple patterns to make the orbital reader robust
ORBITAL_PATTERNS = [
    r"THIS SYSTEM HAS\s+(\d+)\s+VALENCE MOLECULAR ORBITALS,\s+(\d+)\s+OF THEM OCCUPIED",
    r"HAS\s+(\d+)\s+VALENCE\s+(?:MOLECULAR\s+)?ORBITALS.*?(\d+)\s+OCCUPIED",
    r"(\d+)\s+VALENCE\s+MO(?:LECULAR)?S?.*?(\d+)\s+OCCUPIED",
    r"VALENCE\s+ORBITALS.*?(\d+).*?(\d+)\s+OCCUPIED"
]


def read_valence_orbitals(lines: List[str]) -> Tuple[Optional[int], Optional[int]]:
    """
    Reads the number of valence molecular orbitals and how many are occupied.
    Supports multiple .ndo output formats.
    """
    for line in lines:
        for pattern in ORBITAL_PATTERNS:
            match = re.search(pattern, line, re.IGNORECASE)
            if match:
                total, occupied = match.groups()
                return int(total), int(occupied)

    return None, None