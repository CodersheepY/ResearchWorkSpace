from ase.io import read, write

# read your CIF file
structure = read('Ba8Zr8O24.cif')

# write the POSCAR file
write('POSCAR1', structure, format='vasp')
