import re
from typing import List, Dict, Any

def leer_estructura_molecular(lines: List[str]) -> List[Dict[str, Any]]:
    estructura = []
    dentro_coords = False
    start_idx = None
    for i, line in enumerate(lines):
        if 'INPUT ATOMIC CARTESIAN COORDINATES (A)' in line:
            dentro_coords = True
            continue
        if dentro_coords and 'ATOM      X         Y         Z      NAT' in line:
            start_idx = i + 2
            break
    if start_idx is None:
        raise ValueError("No se encontró bloque de coordenadas.")
    for line in lines[start_idx:]:
        if line.strip() == '' or not re.match(r'\s*[A-Z]', line):
            break
        parts = line.split()
        if len(parts) >= 6:
            simbolo = parts[0]
            numero = int(parts[1])
            x, y, z = map(float, parts[2:5])
            nat = int(parts[5])
            valencia = {'H': 1, 'C': 4, 'O': 4}.get(simbolo, 0)
            estructura.append({
                'simbolo': simbolo,
                'numero': numero,
                'x': x,
                'y': y,
                'z': z,
                'nat': nat,
                'valencia': valencia
            })
    return estructura