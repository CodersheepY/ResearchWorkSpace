from ase.build import bulk
from ase.calculators.vasp import Vasp
from ase.dft.bandgap import bandgap

si = bulk(name="Si", crystalstructure="diamond", a=5.43)
si = si*[2, 2, 2]
si.calc = Vasp(prec="normal", xc="pbe",
                encut=400.0, kpts=[2, 2, 2], ismear=0)

energy = si.get_potential_energy()
print(f"Energy = {energy:5.3f} eV")

gap, p1, p2 = bandgap(si.calc, direct=True)
print(f"gap = {gap:5.3f}, p1 = {p1}, p2 = {p2}")