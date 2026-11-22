"""Configuración de corridas no interactivas."""
from dataclasses import dataclass
from pathlib import Path

from .ndol_parser import EstadoCI, NDOLOutput


@dataclass
class Configuracion:
    ndol: Path
    output_dir: Path
    export_ao: bool
    estados: dict[str, list[int] | str]
    fragmentos: dict[str, str]


def leer_configuracion(ruta: Path) -> Configuracion:
    try:
        import tomllib
    except ImportError:
        try:
            import tomli as tomllib
        except ImportError as exc:
            raise ValueError("TOML on Python 3.10 requires tomli. Install it with python -m pip install tomli.") from exc

    ruta = ruta.resolve()
    with ruta.open("rb") as archivo:
        datos = tomllib.load(archivo)

    desconocidas = datos.keys() - {"ndol", "output_dir", "export_ao", "states", "fragments"}
    if desconocidas:
        raise ValueError("Unknown configuration keys: " + ", ".join(sorted(desconocidas)))

    rutas = {}
    for clave in ("ndol", "output_dir"):
        valor = datos.get(clave)
        if not isinstance(valor, str) or not valor.strip():
            raise ValueError(f"{clave} must be a nonempty path string.")
        rutas[clave] = (ruta.parent / valor).resolve()

    export_ao = datos.get("export_ao", False)
    if not isinstance(export_ao, bool):
        raise ValueError("export_ao must be true or false.")

    estados = datos.get("states")
    if not isinstance(estados, dict) or not estados or estados.keys() - {"singlet", "triplet"}:
        raise ValueError("[states] must specify singlet, triplet or both.")
    for espin, seleccion in estados.items():
        if seleccion == "all":
            continue
        if (not isinstance(seleccion, list) or not seleccion
                or any(type(n) is not int or n < 1 for n in seleccion)
                or len(set(seleccion)) != len(seleccion)):
            raise ValueError(f'states.{espin} must be "all" or a nonempty list of distinct positive state numbers.')

    fragmentos = datos.get("fragments")
    if (not isinstance(fragmentos, dict) or not fragmentos
            or any(not nombre.strip() or not isinstance(valor, str) or not valor.strip()
                   for nombre, valor in fragmentos.items())):
        raise ValueError('[fragments] must map fragment names to ATOM number strings, such as "1-4, 8".')

    return Configuracion(rutas["ndol"], rutas["output_dir"], export_ao, estados, fragmentos)


def seleccionar_estados(ndol: NDOLOutput, config: Configuracion) -> list[EstadoCI]:
    encontrados = []
    for espin, disponibles in (
        ("singlet", ndol.estados_singlete), ("triplet", ndol.estados_triplete)
    ):
        seleccion = config.estados.get(espin)
        if seleccion is None:
            continue
        if not disponibles:
            raise ValueError(f"No {espin} states were found in the NDOL output.")
        if seleccion == "all":
            encontrados.extend(disponibles)
            continue
        por_numero = {estado.indice: estado for estado in disponibles}
        faltantes = [n for n in seleccion if n not in por_numero]
        if faltantes:
            raise ValueError(f"Unknown {espin} NDOL state numbers: " + ", ".join(map(str, faltantes)))
        encontrados.extend(por_numero[n] for n in seleccion)
    return encontrados
