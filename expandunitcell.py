from pymatgen.core import Structure

# Load the BaAlO₃ CIF file
structure = Structure.from_file("BaFeO3.cif")

# Expand the unit cell by 2x2x2 to get Ba8Al8O24
expanded_structure = structure * (2, 2, 2)

# Save the expanded structure as a new CIF file
expanded_structure.to(filename="Ba8Fe8O24.cif")
