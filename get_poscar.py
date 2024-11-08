from ase.io import read, write

# read your CIF file
structure = read('Ba8Fe8O24.cif')

# write the POSCAR file
write('POSCAR2', structure, format='vasp')
