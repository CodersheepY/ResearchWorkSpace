from ase.io import read, write

# read your CIF file
structure = read('BaZrO3.cif')

# write the POSCAR file
write('POSCAR4', structure, format='vasp')
