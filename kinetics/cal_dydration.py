def setup_directory(structure_name):
    """
    Creates a directory if it does not exist.
    """
    import os
    if not os.path.exists(structure_name):
        os.makedirs(structure_name)

def get_hydration_energy(cif_file):
    """
    Calculates and prints defect formation and hydration energies with ZPE corrections.

    Args:
        cif_file (str): Path to the CIF file.
    Returns:
        delta_E_hydr (float): Hydration energy.
    """

    from ase.io import read
    from ase import Atom
    from ase.io.vasp import read_vasp_out

    # --- Constant energy values for H2O, H2, O2 (eV)
    E_H2O = -14.891
    # E_H2 = -6.715
    # E_O2 = -9.896

    # --- Constant ZPE values (eV)
    ZPE_H2O = 0.56
    # ZPE_H2 = 0.27
    # ZPE_O2 = 0.10

    # --- Apply ZPE corrections
    # E_O2_corr = 2 * ((E_H2O + ZPE_H2O) - (E_H2 + ZPE_H2)) - ZPE_O2
    # E_H2_corr = E_H2 + ZPE_H2
    E_H2O_corr = E_H2O + ZPE_H2O
    # energy_of_H2O_molecule = E_H2O + ZPE_H2O

    # --- Read the CIF file
    atoms = read(cif_file)

    # --- make the hydrated structure (hydrated)
    # --- Copy the original structure
    atoms_with_Vox = atoms.copy()

    # --- Step1: Remove an oxygen atom to simulate the vacancy(Vox)
    # --- Find the oxygen atom 
    oxygen_indices = [atom.index for atom in atoms_with_Vox if atom.symbol == "O"]
    # --- Record the oxygen atom position
    oxygen_position = atoms[oxygen_indices[0]].position
    # --- Remove an oxygen atom
    minus_O = atoms_with_Vox.pop(oxygen_indices[0])

    # --- Step2: Add back the oxygen atom and a hydrogen atom to simulate proton defect(OHx)
    # --- Copy the structure with the vacancy
    atoms_with_OHx=atoms_with_Vox.copy()
    # --- Add an oxygen atom
    atoms_with_OHx.append(Atom("O", position=oxygen_position))
    # --- Add a hydrogen atom
    atoms_with_OHx.append(Atom("H", position=oxygen_position + [0.0, 0.0, 0.96]))

    # --- Set up directories for VASP runs
    setup_directory("pristine")
    setup_directory("Vox")
    setup_directory("OHx")

    # --- Write atoms to corresponding directories
    atoms.write("pristine/POSCAR")
    atoms_with_Vox.write("Vox/POSCAR")
    atoms_with_OHx.write("OHx/POSCAR")
    
    # --- Ensure the OUTCAR files exist
    for dir_name in ['pristine', 'Vox', 'OHx']:
        if not os.path.exists(f'{dir_name}/OUTCAR'):
            raise FileNotFoundError(f"OUTCAR not found in {dir_name} directory. Make sure VASP has finished running.")
    
    # --- Retrieve the energies from the OUTCAR files
    atoms_pristine = read_vasp_out('pristine/OUTCAR')
    E_pristine = atoms_pristine.get_potential_energy()

    atoms_Vox = read_vasp_out('Vox/OUTCAR')
    E_Vox = atoms_Vox.get_potential_energy()

    atoms_OHx = read_vasp_out('OHx/OUTCAR')
    E_OHx = atoms_OHx.get_potential_energy()

    # --- Calculate hydration energy
    # delta_E_Vox = E_Vox - E_pristine + E_O2_corr / 2
    # delta_E_OHx = E_OHx - E_pristine - ((E_H2O_corr - E_H2_corr / 2) / 2)
    # delta_E_hydr = (2*energy_of_plus_H) - (energy_of_original + energy_of_minus_O + energy_of_H2O_molecule)
    delta_E_hydr = (2 * E_OHx) - (E_pristine + E_Vox + E_H2O_corr)

    print(f"E_pristine: {E_pristine}, E_Vox: {E_Vox}, E_OHx: {E_OHx}")
    print(f"Calculated hydration energy: {delta_E_hydr} eV")

    # Write hydration energy to vasp.out
    with open('vasp.out', 'a') as f:
        f.write(f'Pristine Energy: {E_pristine} eV\n')
        f.write(f'Vox Energy: {E_Vox} eV\n')
        f.write(f'OHx Energy: {E_OHx} eV\n')
        f.write(f'Hydration Energy: {delta_E_hydr} eV\n')

    # --- Return results
    return delta_E_hydr


if __name__ == "__main__":
    get_hydration_energy("BaZrO3.cif")
    