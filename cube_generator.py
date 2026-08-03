import numpy as np
from typing import List, Dict, Tuple, Optional

class CubeGenerator:
    """
    Módulo para generar mapas 3D de densidad de transición, hueco y partícula en formato .cube
    """
    
    def __init__(self):
        self.grid_info = None
        
    def build_transition_matrices(self, cis_data: Dict, n_occ: int, n_virt: int, state_index: int = 0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Construye la matriz de transición T y las densidades de hueco/partícula en base MO.
        
        Args:
            cis_data: Datos CIS parseados
            n_occ: Número de orbitales ocupados
            n_virt: Número de orbitales virtuales
            state_index: Índice del estado excitado
            
        Returns:
            Tuple con (T, D_h_MO, D_p_MO)
        """
        states = cis_data["excited_states"]
        state = states[state_index]
        T = np.zeros((n_occ, n_virt))
        
        print(f"  Construyendo matriz T para estado {state_index + 1}")
        print(f"  Dimensiones: {n_occ} ocupados × {n_virt} virtuales")
        
        transitions_used = 0
        transitions_skipped = 0
        
        # ANALIZAR los índices de orbitales para entender el problema
        occupied_indices = []
        virtual_indices = []
        
        for tr in state["transitions"]:
            occupied_indices.append(tr["occupied_orb"])
            virtual_indices.append(tr["virtual_orb"])
        
        if occupied_indices:
            print(f"  Índices ocupados encontrados: {min(occupied_indices)} - {max(occupied_indices)}")
            print(f"  Índices virtuales encontrados: {min(virtual_indices)} - {max(virtual_indices)}")
        
        for tr in state["transitions"]:
            i = tr["occupied_orb"] - 1  # Convertir a índice 0-based
            a = tr["virtual_orb"] - 1   # Convertir a índice 0-based
            
            # VERIFICAR LÍMITES antes de asignar
            if i < 0 or i >= n_occ:
                print(f"    Advertencia: índice ocupado {i+1} fuera de rango [1, {n_occ}]")
                transitions_skipped += 1
                continue
            if a < 0 or a >= n_virt:
                print(f"    Advertencia: índice virtual {a+1} fuera de rango [1, {n_virt}]")
                transitions_skipped += 1
                continue
                
            T[i, a] = tr["coef"]
            transitions_used += 1
        
        print(f"  Transiciones: {transitions_used} usadas, {transitions_skipped} omitidas")
        
        # Solo proceder si tenemos transiciones válidas
        if transitions_used > 0:
            D_h_MO = T @ T.T  # Densidad de hueco en base MO
            D_p_MO = T.T @ T  # Densidad de partícula en base MO
            print(f"  ✓ Matrices de densidad construidas")
        else:
            print(f"  ✗ ERROR: No hay transiciones válidas para el estado {state_index + 1}")
            print(f"    Los índices de orbitales en CIS no coinciden con las dimensiones MO")
            print(f"    Se necesitan ajustes en la numeración de orbitales")
            # Crear matrices cero
            D_h_MO = np.zeros((n_occ, n_occ))
            D_p_MO = np.zeros((n_virt, n_virt))
        
        return T, D_h_MO, D_p_MO

    def build_transition_matrices_with_offset(self, cis_data: Dict, n_occ: int, n_virt: int, 
                                           occupied_offset: int, virtual_offset: int, 
                                           state_index: int = 0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Versión alternativa que usa offsets para corregir la numeración de orbitales.
        
        Args:
            cis_data: Datos CIS parseados
            n_occ: Número de orbitales ocupados
            n_virt: Número de orbitales virtuales
            occupied_offset: Offset para orbitales ocupados (primer orbital ocupado - 1)
            virtual_offset: Offset para orbitales virtuales (primer orbital virtual - 1)
            state_index: Índice del estado excitado
        """
        states = cis_data["excited_states"]
        state = states[state_index]
        T = np.zeros((n_occ, n_virt))
        
        print(f"  Construyendo matriz T con offsets para estado {state_index + 1}")
        print(f"  Offsets: ocupados={occupied_offset}, virtuales={virtual_offset}")
        
        transitions_used = 0
        
        for tr in state["transitions"]:
            i = tr["occupied_orb"] - occupied_offset - 1  # Aplicar offset y convertir a 0-based
            a = tr["virtual_orb"] - virtual_offset - 1    # Aplicar offset y convertir a 0-based
            
            # VERIFICAR LÍMITES después de aplicar offset
            if i < 0 or i >= n_occ:
                continue
            if a < 0 or a >= n_virt:
                continue
                
            T[i, a] = tr["coef"]
            transitions_used += 1
        
        print(f"  Transiciones usadas con offset: {transitions_used}")
        
        if transitions_used > 0:
            D_h_MO = T @ T.T
            D_p_MO = T.T @ T
            print(f"  ✓ Matrices de densidad construidas con offsets")
            return T, D_h_MO, D_p_MO
        else:
            print(f"  ✗ ERROR: No hay transiciones válidas incluso con offsets")
            return T, np.zeros((n_occ, n_occ)), np.zeros((n_virt, n_virt))

    def transform_to_AO(self, D_MO: np.ndarray, C_MO: np.ndarray) -> np.ndarray:
        """
        Proyecta matriz de densidad desde base MO a base AO.
        """
        return C_MO @ D_MO @ C_MO.T

    def create_grid_around_molecule(self, structure: List[Dict], padding: float = 5.0, 
                                   points_per_angstrom: int = 10) -> Dict:
        """
        Crea una grilla 3D que engloba toda la molécula.
        """
        # Extraer coordenadas
        coords = np.array([[atom["x"], atom["y"], atom["z"]] for atom in structure])
        
        # Encontrar límites
        min_coords = coords.min(axis=0) - padding
        max_coords = coords.max(axis=0) + padding
        
        # Calcular número de puntos
        ranges = max_coords - min_coords
        nx = int(ranges[0] * points_per_angstrom)
        ny = int(ranges[1] * points_per_angstrom)
        nz = int(ranges[2] * points_per_angstrom)
        
        # Asegurar un mínimo de puntos
        nx = max(nx, 10)
        ny = max(ny, 10)
        nz = max(nz, 10)
        
        # Espaciamiento
        dx = ranges[0] / (nx - 1) if nx > 1 else 1.0
        dy = ranges[1] / (ny - 1) if ny > 1 else 1.0
        dz = ranges[2] / (nz - 1) if nz > 1 else 1.0
        
        # Generar puntos de la grilla
        grid_points = []
        for i in range(nx):
            x = min_coords[0] + i * dx
            for j in range(ny):
                y = min_coords[1] + j * dy
                for k in range(nz):
                    z = min_coords[2] + k * dz
                    grid_points.append((x, y, z))
        
        grid_info = {
            "points": grid_points,
            "shape": (nx, ny, nz),
            "origin": min_coords,
            "spacing": (dx, dy, dz),
            "min_coords": min_coords,
            "max_coords": max_coords
        }
        
        self.grid_info = grid_info
        return grid_info

    def evaluate_density_on_grid(self, D_AO: np.ndarray, basis_functions: List, 
                               grid_points: List[Tuple]) -> np.ndarray:
        """
        Evalúa densidad en cada punto del grid.
        """
        values = np.zeros(len(grid_points))
        
        print(f"    Evaluando en {len(grid_points)} puntos...")
        
        for idx, r in enumerate(grid_points):
            if idx % 10000 == 0:  # Progress indicator
                print(f"      Progreso: {idx}/{len(grid_points)} puntos")
                
            density_val = 0.0
            r_array = np.array(r)
            
            # Evaluar todas las funciones de base en este punto
            basis_vals = np.array([bf(r_array) for bf in basis_functions])
            
            # Calcular densidad: sum_{mu,nu} D_AO[mu,nu] * chi_mu(r) * chi_nu(r)
            density_val = np.einsum('i,ij,j->', basis_vals, D_AO, basis_vals)
            
            values[idx] = density_val
            
        return values

    def write_cube_file(self, filename: str, grid_shape: Tuple, origin: Tuple, 
                       spacing: Tuple, values: np.ndarray, structure: List[Dict]):
        """
        Escribe archivo .cube con densidad evaluada.
        """
        nx, ny, nz = grid_shape
        
        with open(filename, 'w') as f:
            # Cabecera
            f.write("CUBE file generated by MONpy\n")
            f.write("Electron/hole density map\n")
            
            # Información de la grilla
            f.write(f"{len(structure):5d}{origin[0]:12.6f}{origin[1]:12.6f}{origin[2]:12.6f}\n")
            f.write(f"{nx:5d}{spacing[0]:12.6f} 0.000000 0.000000\n")
            f.write(f"{ny:5d} 0.000000{spacing[1]:12.6f} 0.000000\n")
            f.write(f"{nz:5d} 0.000000 0.000000{spacing[2]:12.6f}\n")
            
            # Información atómica
            for atom in structure:
                f.write(f"{atom['number']:5d} 0.000000 {atom['x']:12.6f} {atom['y']:12.6f} {atom['z']:12.6f}\n")
            
            # Datos de densidad
            vals_3d = values.reshape(grid_shape)
            flat_vals = vals_3d.flatten()
            
            # Escribir en bloques de 6 valores por línea
            for i in range(0, len(flat_vals), 6):
                line_vals = flat_vals[i:i+6]
                f.write(" ".join(f"{v:13.5E}" for v in line_vals) + "\n")

    def generate_density_maps(self, structure: List[Dict], cis_data: Dict, 
                            C_occ: np.ndarray, C_virt: np.ndarray, 
                            basis_functions: List, state_index: int = 0,
                            padding: float = 5.0, resolution: int = 10,
                            use_offset: bool = False, occupied_offset: int = 300, 
                            virtual_offset: int = 300) -> Dict:
        """
        Función principal para generar todos los mapas de densidad.
        
        Args:
            use_offset: Si True, usa offsets para corregir numeración de orbitales
            occupied_offset: Offset para orbitales ocupados
            virtual_offset: Offset para orbitales virtuales
        """
        print(f"  Generando mapas 3D para estado excitado {state_index + 1}...")
        
        # Determinar dimensiones
        n_occ = C_occ.shape[1]
        n_virt = C_virt.shape[1]
        
        print(f"  Dimensiones MO: {n_occ} ocupados, {n_virt} virtuales")
        print(f"  Dimensiones AO: {C_occ.shape[0]} orbitales atómicos")
        
        # Construir matrices de transición
        if use_offset:
            T, D_h_MO, D_p_MO = self.build_transition_matrices_with_offset(
                cis_data, n_occ, n_virt, occupied_offset, virtual_offset, state_index)
        else:
            T, D_h_MO, D_p_MO = self.build_transition_matrices(cis_data, n_occ, n_virt, state_index)
        
        # Verificar si hay transiciones válidas
        if np.all(T == 0):
            print(f"  ✗ No se pueden generar mapas para estado {state_index + 1}: sin transiciones válidas")
            return None
        
        # Proyectar a base AO
        print("  Proyectando a base AO...")
        D_h_AO = self.transform_to_AO(D_h_MO, C_occ)
        D_p_AO = self.transform_to_AO(D_p_MO, C_virt)
        
        # Densidad de transición total
        D_trans_AO = D_h_AO + D_p_AO
        
        # Crear grilla
        print("  Creando grilla 3D...")
        grid_info = self.create_grid_around_molecule(structure, padding, resolution)
        
        # Evaluar densidades en la grilla
        print("  Evaluando densidades en la grilla...")
        rho_h = self.evaluate_density_on_grid(D_h_AO, basis_functions, grid_info["points"])
        rho_p = self.evaluate_density_on_grid(D_p_AO, basis_functions, grid_info["points"])
        rho_trans = self.evaluate_density_on_grid(D_trans_AO, basis_functions, grid_info["points"])
        
        # Escribir archivos .cube
        base_name = f"state_{state_index + 1}"
        
        self.write_cube_file(f"{base_name}_hole.cube", grid_info["shape"], 
                           grid_info["origin"], grid_info["spacing"], rho_h, structure)
        self.write_cube_file(f"{base_name}_particle.cube", grid_info["shape"], 
                           grid_info["origin"], grid_info["spacing"], rho_p, structure)
        self.write_cube_file(f"{base_name}_transition.cube", grid_info["shape"], 
                           grid_info["origin"], grid_info["spacing"], rho_trans, structure)
        
        print(f"  ✓ Mapas 3D generados exitosamente:")
        print(f"     - {base_name}_hole.cube (densidad de hueco)")
        print(f"     - {base_name}_particle.cube (densidad de partícula)")
        print(f"     - {base_name}_transition.cube (densidad de transición total)")
        
        return {
            "hole_density": f"{base_name}_hole.cube",
            "particle_density": f"{base_name}_particle.cube",
            "transition_density": f"{base_name}_transition.cube",
            "grid_info": grid_info
        }

# Funciones auxiliares para bases atómicas (sin cambios)
def create_s_orbital(center: np.ndarray, alpha: float):
    def s_orbital(r):
        r_diff = r - center
        return np.exp(-alpha * np.dot(r_diff, r_diff))
    return s_orbital

def create_p_orbital(center: np.ndarray, alpha: float, axis: int):
    def p_orbital(r):
        r_diff = r - center
        r_sq = np.dot(r_diff, r_diff)
        return r_diff[axis] * np.exp(-alpha * r_sq)
    return p_orbital

def create_minimal_basis(structure: List[Dict]) -> List:
    basis_functions = []
    basis_sets = {
        1: [('s', 0.5)],   # H
        6: [('s', 0.5), ('p', 0.5), ('p', 0.5), ('p', 0.5)],  # C
        7: [('s', 0.5), ('p', 0.5), ('p', 0.5), ('p', 0.5)],  # N
        8: [('s', 0.5), ('p', 0.5), ('p', 0.5), ('p', 0.5)],  # O
    }
    
    for atom in structure:
        Z = atom["number"]
        center = np.array([atom["x"], atom["y"], atom["z"]])
        
        if Z in basis_sets:
            for orb_type, alpha in basis_sets[Z]:
                if orb_type == 's':
                    basis_functions.append(create_s_orbital(center, alpha))
                elif orb_type == 'p':
                    for axis in range(3):
                        basis_functions.append(create_p_orbital(center, alpha, axis))
        else:
            basis_functions.append(create_s_orbital(center, 0.5))
    
    print(f"  Base mínima creada: {len(basis_functions)} funciones")
    return basis_functions