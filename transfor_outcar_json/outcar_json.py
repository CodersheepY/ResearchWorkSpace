from ase.io.vasp import read_vasp_out
import json
from pymatgen.core import Composition

def extract_vasp_data_with_fixed_formula(outcar_path, output_json_path):
    """
    Extract energy, forces, and stress from the OUTCAR file and save to a JSON file with a fixed chemical formula order.

    Args:
    - outcar_path: str, path to the OUTCAR file
    - output_json_path: str, path to the output JSON file
    """
    try:
        # Read the OUTCAR file
        atoms = read_vasp_out(outcar_path)

        # Extract and standardize the chemical formula
        raw_formula = atoms.get_chemical_formula(empirical=True)
        standardized_formula = Composition(raw_formula).reduced_formula  # Standardized chemical formula

        # Extract energy, forces, and stress
        energy = atoms.get_potential_energy()  # eV
        forces = atoms.get_forces()           # eV/Å
        stress = atoms.get_stress()           # kBar, convert to GPa
        stress_gpa = [s / 10 for s in stress]

        # Construct the data dictionary
        data = {
            "material": standardized_formula,  
            "energy (eV)": energy,
            "forces (eV/Å)": forces.tolist(),
            "stress (GPa)": stress_gpa
        }

        # Save to the JSON file
        with open(output_json_path, 'w') as f:
            json.dump(data, f, indent=4)

        print(f"Data for {standardized_formula} successfully saved to {output_json_path}!")

    except Exception as e:
        print(f"An error occurred: {e}")

# Set the input and output file paths
outcar_path = "OUTCAR"          
output_json_path = "output.json" 

# Extract data and save
extract_vasp_data_with_fixed_formula(outcar_path, output_json_path)