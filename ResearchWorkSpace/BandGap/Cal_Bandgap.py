from ase import Atoms
from ase.calculators.espresso import Espresso
from ase.dft.bandgap import bandgap

# Create a silicon crystal
silicon = Atoms('Si', positions=[(0, 0, 0)], cell=[5.43, 5.43, 5.43], pbc=True)

# Create a calculator object
calc = Espresso(pseudopotentials={'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'},
                input_data={
                    'control': {'calculation': 'scf'},  # 自洽场计算
                    'system': {'ecutwfc': 30},  # 截断能量（简化设置）
                })

# Attach the calculator to the silicon crystal
silicon.calc = calc

# Calculate the band gap
gap, p1, p2 = bandgap(silicon.calc)
print(f"Band gap: {gap} eV")
