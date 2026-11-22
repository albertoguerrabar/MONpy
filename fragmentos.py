import re

from ndol_parser import Atomo


def definir_fragmentos(atomos: list[Atomo]) -> dict[str, list[int]]:
    numeros = {atomo.numero for atomo in atomos}
    print("\nATOM numbers: " + ", ".join(map(str, sorted(numeros))))
    while True:
        cantidad = input("Number of fragments: ").strip()
        if cantidad.isdecimal() and 1 <= int(cantidad) <= len(numeros):
            break
        print(f"Enter a number between 1 and {len(numeros)}.")

    fragmentos: dict[str, list[int]] = {}
    disponibles = numeros.copy()
    total = int(cantidad)
    for i in range(total):
        while True:
            nombre = input(f"Fragment name {i + 1}: ").strip()
            if nombre and nombre not in fragmentos:
                break
            print("Enter a nonempty, unique fragment name.")

        while True:
            entrada = input(f"ATOM numbers for '{nombre}' (e.g. 1,3-6): ").strip()
            seleccion = []
            for campo in entrada.split(","):
                intervalo = re.fullmatch(r"(\d+)(?:\s*-\s*(\d+))?", campo.strip())
                if intervalo is None:
                    break
                inicio = int(intervalo[1])
                fin = int(intervalo[2] or inicio)
                if not 0 <= fin - inicio < len(numeros):
                    break
                seleccion.extend(range(inicio, fin + 1))
            else:
                elegidos = set(seleccion)
                if len(elegidos) != len(seleccion) or not elegidos <= disponibles:
                    print("Use existing ATOM numbers once each.")
                    continue
                restantes = disponibles - elegidos
                pendientes = total - i - 1
                if len(restantes) < pendientes or (not pendientes and restantes):
                    print("Assign every atom and leave at least one for each remaining fragment.")
                    print("Unassigned after this selection: " + ", ".join(map(str, sorted(restantes))))
                    continue
                fragmentos[nombre] = seleccion
                disponibles = restantes
                break
            print("Enter ATOM numbers or ascending ranges separated by commas.")
    return fragmentos
