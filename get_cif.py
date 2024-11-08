from pymatgen.ext.matproj import MPRester

# Initialize the MP Rester
mpr = MPRester("kzum4sPsW7GCRwtOqgDIr3zhYrfpaguK")
structure = mpr.get_structure_by_material_id("mp-18965") # Input the material ID
structure.to(filename="BaCoO3.cif")
