from ase.io import read

def calculate_oh_bond_lengths(poscar_path: str, oh_tolerance: float = 1.2):
    """
    Read the POSCAR file and calculate all O-H bond lengths.

    Args:
        poscar_path (str): Path to the POSCAR file.
        oh_tolerance (float): Maximum O-H bond length threshold (Å).

    Returns:
        list: A list containing each O-H bond and its distance.
    """
    # Read the structure
    atoms = read(poscar_path, format='vasp')
    
    # Get the indices of O and H atoms
    o_indices = [i for i, symbol in enumerate(atoms.get_chemical_symbols()) if symbol == 'O']
    h_indices = [i for i, symbol in enumerate(atoms.get_chemical_symbols()) if symbol == 'H']
    
    oh_bonds = []
    
    # Iterate over all O-H combinations and calculate distances
    for o_idx in o_indices:
        for h_idx in h_indices:
            dist = atoms.get_distance(o_idx, h_idx, mic=True)  # mic=True considers periodic boundary conditions
            if dist <= oh_tolerance:
                oh_bonds.append((o_idx, h_idx, dist))
    
    return oh_bonds

# Use the function
poscar_file = "Ba8Zr8O24_H1.vasp"  # Replace with the actual POSCAR file path
oh_bond_lengths = calculate_oh_bond_lengths(poscar_file)

# Output the results
print("O-H bond length analysis:")
for o_idx, h_idx, dist in oh_bond_lengths:
    print(f"O atom index: {o_idx}, H atom index: {h_idx}, Bond length: {dist:.3f} Å")