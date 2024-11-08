from pymatgen.ext.matproj import MPRester

# Initialize the MP Rester
mpr = MPRester("kzum4sPsW7GCRwtOqgDIr3zhYrfpaguK")
structure = mpr.get_structure_by_material_id("mp-1192651") # Input the material ID
structure.to(filename="BaFeO3.cif")
