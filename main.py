# main.py
from parser import read_molecular_structure, read_valence_orbitals
from cis import read_cis_transitions
from densidad_transicion import build_fragment_transition_density
from propiedades import calculate_all_properties_for_state
from graficar_densidad import plot_transition_density
from fragmentos import define_fragments
from cube_generator import CubeGenerator, create_minimal_basis

import numpy as np
import os

LOGO = """
           __  __  ____  _   _             
          |  \/  |/ __ \| \ | |            
          | \  / | |  | |  \| |_ __  _   _ 
          | |\/| | |  | | . ` | '_ \| | | |
          | |  | | |__| | |\  | |_} | |_| |
          |_|  |_|\____/|_| \_| .__/ \__, |
                              | |     __/ |
                              |_|    |___/ 
"""

INFO_BLOCK = """
 *******************************************************
 *               Alberto Guerra-Barroso                *
 *    Lab of Computational and Theoretical Chemistry   *
 *      Faculty of Chemistry,University of Havana      *
 *                                                     *
 *                        and                          *
 *                                                     *
 *              Ana Lilian Montero-Alejo               *
 *          Department of Physics, FCNMM, UTEM         *                    
 *******************************************************
"""

def write_results(filename: str, logo: str, info_block: str, structure, fragments, results):
    with open(filename, "w") as f:
        f.write(logo + "\n")
        f.write(info_block + "\n\n")
        f.write("=== MONpy analysis results ===\n\n")
        f.write(f"Number of atoms: {len(structure)}\n")
        f.write(f"Number of fragments: {len(fragments)}\n\n")
        f.write("Fragments:\n")
        for fid, atoms in fragments.items():
            f.write(f"  Fragment {fid}: {atoms}\n")
        f.write("\n")
        for state_idx, res in results.items():
            f.write(f"--- Excited state {state_idx+1} ---\n")
            
            # ✅ ACTUALIZADO: Nuevas propiedades de excitón
            f.write(f"Exciton size (d_exc, Å): {res['properties']['exciton_size_aprox']:.6f}\n")
            f.write(f"Exciton type: {res['properties']['exciton_size_interpretation']}\n")
            
            f.write(f"CTI absolute: {res['properties']['cti_absolute']:.6f}\n")
            f.write(f"CTI Frobenius: {res['properties']['cti_frobenius']:.6f}\n")
            f.write("\nIFCT summary:\n")
            ifct = res['properties']['ifct']
            f.write(f"  Intrinsic CT%: {ifct['intrinsic_CT_percent']:.4f}\n")
            f.write(f"  Intrinsic LE%: {ifct['intrinsic_LE_percent']:.4f}\n")
            if ifct['apparent_CT_percent'] is not None:
                f.write(f"  Apparent CT%: {ifct['apparent_CT_percent']:.4f}\n")
            
            # Información de normalización
            if 'normalization_info' in res:
                f.write("\nNormalization info:\n")
                norm_info = res['normalization_info']
                f.write(f"  Total contributions: {norm_info['total_contributions']:.6f}\n")
                f.write(f"  M sum: {norm_info['M_sum']:.6f}\n")
                f.write(f"  M_norm sum: {norm_info['M_norm_sum']:.6f}\n")
            
            f.write("\n")
    print(f"Results written to {filename}")

def compute_fragment_centroids(structure, fragments):
    # returns array shape (n_frag, 3) with centroid coordinates
    frag_order = list(fragments.keys())
    centroids = []
    for fid in frag_order:
        atom_indices = fragments[fid]
        coords = []
        for at in structure:
            if at["number"] in atom_indices:
                coords.append([at["x"], at["y"], at["z"]])
        if len(coords) == 0:
            centroids.append(np.zeros(3))
        else:
            centroids.append(np.mean(np.array(coords), axis=0))
    return np.array(centroids)

def get_molecular_orbitals_from_file(lines):
    """
    Intenta extraer coeficientes MO del archivo de salida.
    ESTA FUNCIÓN DEBE SER ADAPTADA SEGÚN EL FORMATO DE TU ARCHIVO .ndo
    """
    print("Buscando coeficientes de orbitales moleculares en el archivo...")
    
    # Aquí debes implementar la lógica específica para tu formato de archivo
    # Esto es un placeholder - necesitas adaptarlo
    
    # Ejemplo simplificado: asumir que tenemos alguna información de orbitales
    n_atoms = len([line for line in lines if "ATOM" in line])
    
    # Para pruebas, crear matrices de ejemplo
    n_ao = n_atoms * 4  # Estimación aproximada: 4 orbitales por átomo
    n_occ = 0
    n_virt = 0
    
    # Buscar información de orbitales ocupados/virtuales
    for i, line in enumerate(lines):
        if "OCCUPIED" in line and "ORBITALS" in line:
            try:
                n_occ = int(line.split()[-1])
            except:
                pass
        elif "VIRTUAL" in line and "ORBITALS" in line:
            try:
                n_virt = int(line.split()[-1])
            except:
                pass
    
    # Si no se encontró información, usar valores por defecto
    if n_occ == 0:
        n_occ = n_atoms * 2  # Estimación aproximada
    if n_virt == 0:
        n_virt = n_atoms * 2  # Estimación aproximada
    
    print(f"Estimados: n_ao={n_ao}, n_occ={n_occ}, n_virt={n_virt}")
    
    # Crear matrices de coeficientes de ejemplo (REEMPLAZAR CON DATOS REALES)
    C_occ = np.random.rand(n_ao, n_occ)
    C_virt = np.random.rand(n_ao, n_virt)
    
    # Normalizar las "funciones de base"
    for i in range(n_occ):
        C_occ[:, i] /= np.linalg.norm(C_occ[:, i])
    for i in range(n_virt):
        C_virt[:, i] /= np.linalg.norm(C_virt[:, i])
    
    return C_occ, C_virt, n_ao

def main():
    print(LOGO)
    print(INFO_BLOCK)
    print("\n=== MONpy (fragment-level TDM) ===\n")

    filename = input("Enter .ndo filename: ").strip()
    if not os.path.isfile(filename):
        print("File not found.")
        return

    with open(filename, 'r', encoding='latin1', errors='ignore') as f:
        lines = f.readlines()

    # parse structure and orbitals
    structure = read_molecular_structure(lines)
    total_valence, occupied_valence = read_valence_orbitals(lines)
    print(f"Read {len(structure)} atoms.")
    if total_valence is None:
        print("Warning: could not read valence orbital info.")

    # read CIS
    cis_data = read_cis_transitions(filename)
    excited_states = cis_data.get("excited_states", [])
    if len(excited_states) == 0:
        print("No excited states found in CIS data.")
        return
    print(f"Found {len(excited_states)} excited states.")

    # fragments (mandatory)
    fragments = define_fragments(structure)

    # compute centroids
    centroids = compute_fragment_centroids(structure, fragments)

    # select states to analyze
    entry = input("Enter excited state indices to analyze (e.g. 1,2 or 'all'): ").strip().lower()
    if entry == "all":
        state_indices = list(range(len(excited_states)))
    else:
        try:
            state_indices = [int(x)-1 for x in entry.split(",")]
        except Exception:
            print("Invalid input.")
            return

    results = {}
    # compute per-state TDM and properties
    for sidx in state_indices:
        if sidx < 0 or sidx >= len(excited_states):
            print(f"State index {sidx+1} out of range; skipping.")
            continue
        print(f"Processing excited state {sidx+1}...")
        
        # Obtener ambas matrices
        M, M_norm = build_fragment_transition_density(structure, fragments, cis_data, state_index=sidx)
        
        # Usar M_norm para calcular propiedades (ya que está normalizada)
        props = calculate_all_properties_for_state(M_norm, centroids)
        
        # Guardar ambas matrices e información de normalización
        results[sidx] = {
            "M": M, 
            "M_norm": M_norm, 
            "properties": props,
            "normalization_info": {
                "total_contributions": np.sum(M),
                "M_sum": np.sum(M),
                "M_norm_sum": np.sum(M_norm)
            }
        }

        # save plot for this state - usar M_norm para graficar (valores normalizados)
        frag_names = [str(fid) for fid in fragments.keys()]
        plot_transition_density(M_norm, vmax_global=None, names=frag_names,
                                title=f"State {sidx+1} TDM (fragment)", save=True,
                                filename=f"state_{sidx+1}_tdm.png")

    # write results
    write_results("results.txt", LOGO, INFO_BLOCK, structure, fragments, results)
    
# 8. Opción de generar mapas 3D .cube
    generate_cube = input("\n¿Generar mapas 3D .cube de densidad? (s/n): ").strip().lower()
    if generate_cube == 's':
        cube_gen = CubeGenerator()
        
        print("Generando mapas 3D...")
        
        # Obtener coeficientes MO
        C_occ, C_virt, n_ao = get_molecular_orbitals_from_file(lines, structure, occupied_valence, total_valence)
        
        # Crear base mínima para visualización
        basis_funcs = create_minimal_basis(structure)
        
        # Preguntar si usar offsets para corregir numeración
        use_offset = input("¿Usar offsets para corregir numeración de orbitales? (s/n): ").strip().lower() == 's'
        occupied_offset = 300
        virtual_offset = 300
        
        if use_offset:
            try:
                occupied_offset = int(input(f"Offset para orbitales ocupados (por defecto {occupied_offset}): ") or occupied_offset)
                virtual_offset = int(input(f"Offset para orbitales virtuales (por defecto {virtual_offset}): ") or virtual_offset)
            except:
                print("Usando valores por defecto para offsets")
        
        cube_files_generated = {}
        
        for sidx in state_indices:
            try:
                print(f"Generando mapas 3D para estado {sidx+1}...")
                cube_files = cube_gen.generate_density_maps(
                    structure, cis_data, C_occ, C_virt, basis_funcs, 
                    state_index=sidx, padding=4.0, resolution=8,
                    use_offset=use_offset, occupied_offset=occupied_offset, 
                    virtual_offset=virtual_offset
                )
                
                if cube_files is not None:
                    print(f"✓ Mapas 3D para estado {sidx+1} generados exitosamente")
                    # Agregar información de archivos .cube a los resultados
                    if sidx in results:
                        results[sidx]['cube_files'] = cube_files
                    cube_files_generated[sidx] = cube_files
                else:
                    print(f"✗ No se pudieron generar mapas 3D para estado {sidx+1}")
                    
            except Exception as e:
                print(f"✗ Error generando mapas 3D para estado {sidx+1}: {e}")
    
    print("\n" + "="*50)
    print("ANÁLISIS COMPLETADO")
    print("="*50)
    print("Archivos generados:")
    print("  - results.txt (resultados del análisis)")
    
    for sidx in state_indices:
        print(f"  - state_{sidx+1}_tdm.png (mapa de densidad de transición)")
        if generate_cube == 's' and sidx in results and 'cube_files' in results[sidx]:
            cube_files = results[sidx]['cube_files']
            # VERIFICAR que cube_files no sea None antes de acceder
            if cube_files is not None:
                print(f"    * {cube_files['hole_density']} (densidad de hueco)")
                print(f"    * {cube_files['particle_density']} (densidad de partícula)")
                print(f"    * {cube_files['transition_density']} (densidad de transición)")
    
    if generate_cube == 's' and cube_files_generated:
        print("\nPara visualizar archivos .cube usa:")
        print("  - VMD: vmd -e archivo.cube")
        print("  - PyMOL: load archivo.cube") 
        print("  - GaussView: Abrir directamente el archivo")
    else:
        print("\nNota: No se generaron archivos .cube debido a problemas con la numeración de orbitales.")
        print("      Los índices en CIS (300+) no coinciden con las dimensiones MO (18).")
        print("      Se necesita ajustar la numeración usando offsets.")

if __name__ == "__main__":
    main()