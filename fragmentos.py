# fragmentos.py

def define_fragments(structure):
    """
    Allows the user to define molecular fragments.
    Accepts atom lists separated by commas, spaces, or any mix of both.
    Ensures all atoms are assigned exactly once.
    """

    total_atoms = len(structure)
    atom_indices = [atom["number"] for atom in structure]

    print("\n=== Fragment Definition ===")
    print(f"The molecule contains {total_atoms} atoms.")
    print("All atoms must belong to exactly one fragment.")
    print("You may separate atom numbers using commas, spaces, or both.\n")

    # Ask for number of fragments
    while True:
        try:
            n_frag = int(input("Enter the number of fragments: "))
            if n_frag <= 0:
                print("Number must be >= 1.")
                continue
            break
        except ValueError:
            print("Invalid input. Please enter an integer.")

    fragments = {}
    assigned = set()

    # Parser that accepts commas, spaces, and mixed
    def parse_atom_list(text):
        text = text.replace(",", " ")   # treat commas as spaces
        parts = text.split()            # split on whitespace
        return [int(x) for x in parts]

    # Input per fragment
    for fid in range(1, n_frag + 1):
        while True:
            raw = input(f"Enter atoms for fragment {fid} (comma or space separated): ").strip()

            try:
                atoms = parse_atom_list(raw)
            except ValueError:
                print("Invalid input. Use integers only.")
                continue

            if not atoms:
                print("A fragment cannot be empty.")
                continue

            # Check atom existence
            if any(a not in atom_indices for a in atoms):
                print("One or more atom indices do not exist in the structure.")
                continue

            # Duplicates inside the same fragment
            if len(atoms) != len(set(atoms)):
                print("Duplicate atoms detected inside the same fragment.")
                continue

            # Check global duplicates
            overlap = assigned.intersection(atoms)
            if overlap:
                print(f"These atoms are already assigned: {sorted(overlap)}")
                continue

            fragments[fid] = atoms
            assigned.update(atoms)
            break

    # Final check
    if set(atom_indices) != assigned:
        missing = sorted(set(atom_indices) - assigned)
        print("\nERROR: Not all atoms were assigned!")
        print("Missing:", missing)
        raise ValueError("Fragment definition incomplete.")

    print("\nFragments successfully defined:")
    for fid, atoms in fragments.items():
        print(f"  Fragment {fid}: {atoms}")

    return fragments