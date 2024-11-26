import numpy as np
from ase import Atoms, Atom
from ase.io import read, write

def find_neighbors(bulk: Atoms, index: int, max_dist: float = 3.0) -> list:
    """
    Find all neighbors of the specified atom within a given distance range.
    """
    neighbors = []
    positions = bulk.positions
    for i, pos in enumerate(positions):
        if i != index:
            dist = bulk.get_distance(index, i, mic=True)
            if dist < max_dist:
                vec = bulk.get_distance(index, i, vector=True, mic=True)
                neighbors.append({'dist': dist, 'vec': vec / np.linalg.norm(vec)})
    return neighbors

def calculate_direction(neighbors: list, max_neighbors: int = 3) -> np.ndarray:
    """
    Calculate the direction vector based on the list of neighbors.
    """
    direction = np.zeros(3)
    for neighbor in sorted(neighbors, key=lambda x: x['dist'])[:max_neighbors]:
        weight = 1.0 / neighbor['dist']  # Weight is inversely proportional to distance
        direction -= neighbor['vec'] * weight
    norm = np.linalg.norm(direction)
    return direction / norm if norm > 1e-6 else np.array([0, 0, 1])

def add_protons_wisely(bulk: Atoms, n_protons: int = 6, oh_bond_length: float = 0.98, neighbor_max_dist: float = 3.0) -> Atoms:
    """Intelligently add H atoms based on the chemical environment"""
    o_indices = [i for i, symbol in enumerate(bulk.symbols) if symbol == 'O']

    for i in range(n_protons):
        available_o = [idx for idx in o_indices
                      if not any(bulk.get_distance(idx, j, mic=True) < 1.2
                      for j, s in enumerate(bulk.symbols) if s == 'H')]

        if not available_o:
            print(f"Warning: Cannot find suitable O atom for H{i+1}")
            continue

        o_idx = available_o[0]
        o_pos = bulk.positions[o_idx]
        neighbors = find_neighbors(bulk, o_idx, max_dist=neighbor_max_dist)
        direction = calculate_direction(neighbors)

        h_pos = o_pos + direction * oh_bond_length
        bulk.append(Atom('H', position=h_pos))
        print(f"Added H{i+1} near O{o_idx} at {h_pos}")

        oh_dist = bulk.get_distance(-1, o_idx, mic=True)
        print(f"OH bond length: {oh_dist:.3f} Å")
        if abs(oh_dist - oh_bond_length) > 0.1:
            print(f"Warning: OH bond length deviates from expected value for H{i+1}")
    
    return bulk

# Use the improved function and save the result
bulk = read("BaZrO3.cif")
bulk = add_protons_wisely(bulk, n_protons=6)

# Save the CIF file with added H atoms
output_filename = "BaZrO3_with_H.cif"
write(output_filename, bulk)
print(f"New CIF file saved as: {output_filename}")

# Calculate the distance between each H and the nearest O atom
for h_idx, h_symbol in enumerate(bulk.symbols):
    if h_symbol == 'H':
        h_pos = bulk.positions[h_idx]
        closest_o_idx = min(
            [(o_idx, bulk.get_distance(h_idx, o_idx)) for o_idx, o_symbol in enumerate(bulk.symbols) if o_symbol == 'O'],
            key=lambda x: x[1]
        )
        print(f"H{h_idx}: Closest O index = {closest_o_idx[0]}, Distance = {closest_o_idx[1]:.3f} Å")