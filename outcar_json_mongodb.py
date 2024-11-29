from ase.io.vasp import read_vasp_out
from pymongo import MongoClient
from pymatgen.core import Composition

def save_outcar_to_mongodb(outcar_path, mongodb_uri, db_name, collection_name):
    """
    Extract data from a single OUTCAR file and save it to MongoDB.

    Args:
    - outcar_path: str, path to the OUTCAR file
    - mongodb_uri: str, MongoDB connection URI
    - db_name: str, database name
    - collection_name: str, collection name
    """
    try:
        # Read the OUTCAR file
        atoms = read_vasp_out(outcar_path)

        # Extract and standardize the chemical formula
        raw_formula = atoms.get_chemical_formula(empirical=True)
        standardized_formula = Composition(raw_formula).reduced_formula

        # Extract data
        energy = atoms.get_potential_energy()  # Energy (eV)
        forces = atoms.get_forces()           # Forces (eV/Å)
        stress = atoms.get_stress()           # Stress (kBar)
        stress_gpa = [s / 10 for s in stress]  # Convert to GPa

        # Construct the data dictionary
        data = {
            "material": standardized_formula,
            "energy_eV": energy,
            "forces_eV_per_A": forces.tolist(),
            "stress_GPa": stress_gpa,
            "source_path": outcar_path
        }

        # Connect to MongoDB
        client = MongoClient(mongodb_uri)
        db = client[db_name]
        collection = db[collection_name]

        # Insert data into MongoDB
        result = collection.insert_one(data)
        print(f"Data for {standardized_formula} saved to MongoDB with ID: {result.inserted_id}")

    except Exception as e:
        print(f"An error occurred: {e}")

outcar_path = "OUTCAR"               # OUTCAR file path
mongodb_uri = "mongodb://localhost:27017"    # MongoDB connection URI
db_name = "vasp_data"                        # Database name
collection_name = "outcar_results"           # Collection name

save_outcar_to_mongodb(outcar_path, mongodb_uri, db_name, collection_name)