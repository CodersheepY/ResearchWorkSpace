from ase.io import read, write

# read your CIF file
structure = read('Ba8Co8O24.cif')

# write the POSCAR file
write('POSCAR3', structure, format='vasp')
