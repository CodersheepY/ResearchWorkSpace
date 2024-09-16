# input data: total energies of pristine, defective, and hydrated systems
E_pristine = -331.28931146  #  pristine material total energy(eV) 
E_Vox = -318.76206100     # Oxygen vacancy defect system total energy 
E_OHx = -333.11950727     # protonic defect system total energy 
E_H2O = -14.891        # H2O total energy
E_H2 = -6.715          # H2 total energy
E_O2 = -9.896         # O2 total energy

# ZPE values obtained from vibrational analysis (example values)
ZPE_H2O = 0.56       # ZPE for water 
ZPE_H2 = 0.27       # ZPE for H2 
ZPE_O2 = 0.10        # ZPE for O2 

E_O2_corr = 2*((E_H2O + ZPE_H2O) -(E_H2 + ZPE_H2)) - ZPE_O2
print(E_O2_corr)
E_H2_corr = E_H2 + ZPE_H2
print(E_H2_corr)
E_H2O_corr = E_H2O + ZPE_H2O
print(E_H2O_corr)

delta_E_Vox = E_Vox - E_pristine + E_O2_corr/2
print(delta_E_Vox)
delta_E_OHx = E_OHx - E_pristine - ((E_H2O_corr - E_H2_corr / 2) / 2)
print(delta_E_OHx)
delta_E_hydr = 2 * E_OHx - E_pristine - E_Vox - E_H2O_corr
print(delta_E_hydr)



# calculate the Oxygen Vacancy and Protonic Defect Formation Energies with ZPE corrections
# E_Vox_corr = delta_E_Vox + E_pristine - E_O2_corr/2
# E_OHx_corr = delta_E_OHx + E_pristine + (E_H2O_corr - E_H2_corr/2)/2


# Output the results with ZPE corrections
print(f"Corrected Defect Formation Energy for Oxygen Vacancy (Vox): {delta_E_Vox} eV")
print(f"Corrected Defect Formation Energy for Protonic Defect (OHx): {delta_E_OHx} eV")
print(f"Corrected Hydration Energy: {delta_E_hydr} eV")

# Write results to file
# with open('Defect_Hydration_Energies.txt', 'w') as out_file:
#     out_file.write(f"Corrected Defect Formation Energy for Oxygen Vacancy (Vox): {delta_Ef_Vox_corr} eV\n")
#     out_file.write(f"Corrected Defect Formation Energy for Protonic Defect (OHx): {delta_Ef_OHx_corr} eV\n")
