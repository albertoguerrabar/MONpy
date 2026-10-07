import re
from typing import List, Tuple, Optional

def leer_orbitales_valencia(lines: List[str]) -> Tuple[Optional[int], Optional[int]]:
    for line in lines:
        match = re.search(r'THIS SYSTEM HAS\s+(\d+)\s+VALENCE MOLECULAR ORBITALS,\s+(\d+)\s+OF THEM OCCUPIED', line)
        if match:
            total = int(match.group(1))
            ocupados = int(match.group(2))
            return total, ocupados
    return None, None