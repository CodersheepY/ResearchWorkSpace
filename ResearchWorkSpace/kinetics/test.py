from ResearchWorkSpace.kinetics.cal_dydration import get_hydration_energy

# define E_pristine, E_Vox, and E_OHx values
E_pristine = -331.28931146  # pristine material total energy (eV)
E_Vox = -318.76206100       # Oxygen vacancy defect system total energy (eV)
E_OHx = -333.11950727       # protonic defect system total energy (eV)

# calculate defect formation and hydration energies
delta_E_Vox, delta_E_OHx, delta_E_hydr = get_hydration_energy(E_pristine, E_Vox, E_OHx)

# Print corrected defect formation energies
print(f"Corrected Defect Formation Energy for Oxygen Vacancy (Vox): {delta_E_Vox:.6f} eV")
print(f"Corrected Defect Formation Energy for Protonic Defect (OHx): {delta_E_OHx:.6f} eV")
print(f"Corrected Hydration Energy: {delta_E_hydr:.6f} eV")
