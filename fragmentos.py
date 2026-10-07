from typing import List, Dict, Optional, Union


def pedir_entero(mensaje: str, min_val: int = 1, max_val: int = 1000) -> int:
    while True:
        try:
            valor = int(input(mensaje))
            if min_val <= valor <= max_val:
                return valor
            print(f"Valor fuera de rango ({min_val}-{max_val})")
        except ValueError:
            print("Introduce un número entero válido.")


def definir_fragmentos(
    estructura: Union[int, list]
) -> Optional[Dict[str, List[int]]]:
    num_atom = estructura if isinstance(estructura, int) else len(estructura)
    print(f"\nNumber of atoms detected: {num_atom}")

    opcion = input(
        "Calculate transition density matrix by fragments? (y/n): "
    ).strip().lower()
    if opcion != 'y':
        return None

    num_frag = pedir_entero(
        "How many fragments do you want to define?:", 1, num_atom
    )

    fragmentos: Dict[str, List[int]] = {}
    for i in range(num_frag):
        nombre = input(f"Fragment name {i+1}: ").strip()
        atomos = input(
            f"Atom numbers for fragment '{nombre}' (comma separated): "
        )
        atomos_lista = [int(a.strip()) for a in atomos.split(',') if a.strip().isdigit()]
        fragmentos[nombre] = atomos_lista
    return fragmentos