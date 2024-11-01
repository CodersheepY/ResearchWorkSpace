#!/usr/bin/env python

########################################### Headers and Imports #################################################
from pymatgen.core import Composition
from pymatgen.entries.computed_entries import ComputedEntry

from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility


from mp_api.client import MPRester
from pymatgen.ext.matproj import MPRester

##################### Global Initialization Hydrogen and Oxygen Experimental Conditions ########################
def initialize_global_variables():
    """
    Initializes global gas entries and chemical potentials for conditions A, C, and X.
    This function sets up the required parameters for each environment (A, C, and X) by
    defining energies, compositions, and entries for relevant gases.
    """

    global H_Ener_A, O_Ener_A, H_Ener_C, O_Ener_C, O_Ener_X, CO2_Ener_X, CO_Ener_X
    global H2_Comp, O2_Comp, H2O_Comp, CO_Comp, CO2_Comp
    global entriesGases_A, entriesGases_C, entriesGases_X
    global locked_Chem_Potential_A, locked_Chem_Potential_C, locked_Chem_Potential_X

    # Define energies for different environments
    H_Ener_A = -4.024
    O_Ener_A = -8.006
    H_Ener_C = -4.997
    O_Ener_C = -6.166
    O_Ener_X = -6.166
    CO2_Ener_X = -25.556
    CO_Ener_X = -20.232

    # Define compositions for gases
    H2_Comp = Composition("H2")
    O2_Comp = Composition("O2")
    H2O_Comp = Composition("H2O")
    CO_Comp = Composition("CO")
    CO2_Comp = Composition("X")

    # Define entries for each condition
    entriesGases_A = [
        ComputedEntry(H2_Comp, H_Ener_A * H2_Comp.num_atoms),
        ComputedEntry(O2_Comp, O_Ener_A * O2_Comp.num_atoms),
        ComputedEntry(H2O_Comp, H_Ener_A * H2_Comp.num_atoms + O_Ener_A * 0.5 * O2_Comp.num_atoms)
    ]
    entriesGases_C = [
        ComputedEntry(H2_Comp, H_Ener_C * H2_Comp.num_atoms),
        ComputedEntry(O2_Comp, O_Ener_C * O2_Comp.num_atoms),
        ComputedEntry(H2O_Comp, H_Ener_C * H2_Comp.num_atoms + O_Ener_C * 0.5 * O2_Comp.num_atoms)
    ]
    entriesGases_X = [
        ComputedEntry(O2_Comp, O_Ener_X * O2_Comp.num_atoms),
        ComputedEntry(CO_Comp, CO_Ener_X),
        ComputedEntry(CO2_Comp, CO2_Ener_X)
    ]

    # Define chemical potentials
    locked_Chem_Potential_A = {'H2': H_Ener_A * 2, 'O2': O_Ener_A * 2}
    locked_Chem_Potential_C = {'O2': O_Ener_C * 2, 'H2': H_Ener_C * 2}
    locked_Chem_Potential_X = {'X': CO2_Ener_X, 'O2': O_Ener_X * 2}

############################################## Material Entry ###################################################
def prepare_material_entries(api_key, TestMat_Comp, TestMat_Ener):
    """
    Prepare material entries and fetch entries from the Materials Project for a given material.

    Args:
        api_key (str): API key for accessing the Materials Project database.
        TestMat_Comp (str): Composition formula of the material.
        TestMat_Ener (float): Energy of the material to be analyzed.

    Returns:
        tuple: Lists of combined entries for conditions A, C, and X, and the ComputedEntry objects for conditions A and C.
    """

    # Define the material's composition
    TestMat_Comp = Composition(TestMat_Comp)

    # Get the number of oxygen atoms in the material composition
    num_O_atoms = TestMat_Comp.get_atomic_fraction('O') * TestMat_Comp.num_atoms

    # Define computed entries for the material under conditions A and C
    TestMat_entry_A = ComputedEntry(TestMat_Comp, TestMat_Ener - O_Ener_A * num_O_atoms)
    TestMat_entry_C = ComputedEntry(TestMat_Comp, TestMat_Ener - O_Ener_C * num_O_atoms)

    # Initialize compatibility module
    compat = MaterialsProject2020Compatibility()
    
    # Automatically extract elements in TestMat_Comp
    elements = [str(element) for element in TestMat_Comp.elements]
    
    # Initialize MPRester to get material entries from the Materials Project
    with MPRester(api_key) as api:
        entries_MP_Org_AC = api.get_entries_in_chemsys(elements + ['O','H']) 
        entries_MP_Org_X = api.get_entries_in_chemsys(elements + ['O','C'])  

    # Process entries with the compatibility module
    entries_MP_Org_AC = compat.process_entries(entries_MP_Org_AC)
    entries_MP_Org_X = compat.process_entries(entries_MP_Org_X)

    # Combine the fetched entries with the VASP computed entries for Conditions A and C
    entries_VASP_A = [TestMat_entry_A]
    entries_VASP_C = [TestMat_entry_C]

    all_entries_A = entries_MP_Org_AC + entries_VASP_A
    all_entries_C = entries_MP_Org_AC + entries_VASP_C
    entriesTotal_X = entries_MP_Org_X  
    
    return all_entries_A, all_entries_C, entriesTotal_X, TestMat_entry_A, TestMat_entry_C
    
############################################## Condition A #####################################################
def calculate_phase_diagram_condition_A(all_entries_A, entriesGases_A, locked_Chem_Potential_A, TestMat_entry_A):
    """
    Calculate the phase diagram and energy above the hull for Condition A (Hydrogen-rich environment).

    Args:
        all_entries_A (list): List of ComputedEntry objects for Condition A.
        entriesGases_A (list): List of ComputedEntry objects for gases under Condition A.
        locked_Chem_Potential_A (dict): Locked chemical potentials for Condition A.
        TestMat_entry_A (ComputedEntry): ComputedEntry object for the material under Condition A.
    
    Returns:
        tuple: Phase diagram for Condition A, energy per atom, and energy above hull.
    """

    # Define species to eliminate
    eliminate_AC = ['H2', 'O2', 'H2O']
    
    # Filter out H2, O2, and H2O from all_entries_A
    H2_entries_A = [e for e in all_entries_A if e.composition.reduced_formula == 'H2']
    O2_entries_A = [e for e in all_entries_A if e.composition.reduced_formula == 'O2']
    H2O_entries_A = [e for e in all_entries_A if e.composition.reduced_formula == 'H2O']
    
    # Filter out these species from all_entries_A
    all_entries_A = list(filter(lambda e: e.composition.reduced_formula not in eliminate_AC, all_entries_A))
    all_entries_A += entriesGases_A

    # Create phase diagram for Condition A
    pd_A = GrandPotentialPhaseDiagram(all_entries_A, locked_Chem_Potential_A)
    energy_per_atom_A = TestMat_entry_A.energy / TestMat_entry_A.composition.num_atoms
    energy_above_hull_A = TestMat_entry_A.energy / TestMat_entry_A.composition.num_atoms - pd_A.get_hull_energy_per_atom(TestMat_entry_A.composition)

    return pd_A, energy_per_atom_A, energy_above_hull_A

############################################## Condition C #####################################################
def calculate_phase_diagram_condition_C(all_entries_C, entriesGases_C, locked_Chem_Potential_C, TestMat_entry_C):
    """
    Calculate the phase diagram and energy above the hull for Condition C (Oxygen-rich environment).

    Args:
        all_entries_C (list): List of ComputedEntry objects for Condition C.
        entriesGases_C (list): List of ComputedEntry objects for gases under Condition C.
        locked_Chem_Potential_C (dict): Locked chemical potentials for Condition C.
        TestMat_entry_C (ComputedEntry): ComputedEntry object for the material under Condition C.

    Returns:
        tuple: Phase diagram for Condition C, energy per atom, and energy above hull.
    """

    # Define species to eliminate
    eliminate_AC = ['H2', 'O2', 'H2O']
    
    # Filter out H2, O2, and H2O from all_entries_C
    H2_entries_C = [e for e in all_entries_C if e.composition.reduced_formula == 'H2']
    O2_entries_C = [e for e in all_entries_C if e.composition.reduced_formula == 'O2']
    H2O_entries_C = [e for e in all_entries_C if e.composition.reduced_formula == 'H2O']
    
    # Filter out these species from all_entries_C
    all_entries_C = list(filter(lambda e: e.composition.reduced_formula not in eliminate_AC, all_entries_C))
    all_entries_C += entriesGases_C

    # Create phase diagram for Condition C
    pd_C = GrandPotentialPhaseDiagram(all_entries_C, locked_Chem_Potential_C)
    energy_per_atom_C = TestMat_entry_C.energy / TestMat_entry_C.composition.num_atoms
    energy_above_hull_C = TestMat_entry_C.energy / 16 - pd_C.get_hull_energy_per_atom(TestMat_entry_C.composition)

    return pd_C, energy_per_atom_C, energy_above_hull_C

############################################## Condition X #####################################################
def calculate_phase_diagram_condition_X(entriesTotal_X, entriesGases_X, locked_Chem_Potential_X, TestMat_Comp, TestMat_Ener):
    """
    Calculate the phase diagram and energy above the hull for Condition X (CO2-rich environment).

    Args:
        entriesTotal_X (list): List of ComputedEntry objects for Condition X.
        entriesGases_X (list): List of ComputedEntry objects for gases under Condition X.
        locked_Chem_Potential_X (dict): Locked chemical potentials for Condition X.
        TestMat_Comp (str): Composition of the material.
        TestMat_Ener (float): Energy of the material.

    Returns:
        tuple: Phase diagram for Condition X, energy per atom, and energy above hull.
    """
    # Define species to eliminate 
    CO_entries_X = [e for e in entriesTotal_X if e.composition.reduced_formula == 'CO']
    CO2_entries_X = [e for e in entriesTotal_X if e.composition.reduced_formula == 'X']
    O2_entries_X = [e for e in entriesTotal_X if e.composition.reduced_formula == 'O2']

    # Filter out CO, CO2, and O2 from entriesTotal_X
    eliminate_X = ['CO', 'X', 'O2']
    all_entries_X = list(filter(lambda e: e.composition.reduced_formula not in eliminate_X, entriesTotal_X))
    
    # Create a ComputedEntry for the material under Condition X
    TestMat_entry_X = ComputedEntry(TestMat_Comp, TestMat_Ener - O_Ener_X * TestMat_Comp.num_atoms)

    # Add back entries for Condition X
    entries_VASP_X = [TestMat_entry_X]
    all_entries_X += entries_VASP_X + entriesGases_X

    # Create phase diagram for Condition X
    pd_X = GrandPotentialPhaseDiagram(all_entries_X, locked_Chem_Potential_X)
    energy_per_atom_X = TestMat_entry_X.energy / TestMat_entry_X.composition.num_atoms
    energy_above_hull_X = TestMat_entry_X.energy / 16 - pd_X.get_hull_energy_per_atom(TestMat_entry_X.composition)

    return pd_X, energy_per_atom_X, energy_above_hull_X

############################################## Main Function ###################################################
def main():
    api_key = "kzum4sPsW7GCRwtOqgDIr3zhYrfpaguK"
    TestMat_Comp = "Ba8Zr8O24"
    TestMat_Ener = -333.71584216

    # Initialize global gas entries and chemical potentials
    initialize_global_variables()

    # Prepare material entries
    all_entries_A, all_entries_C, entriesTotal_X, TestMat_entry_A, TestMat_entry_C = prepare_material_entries(
        api_key, TestMat_Comp, TestMat_Ener
    )

    # Calculate phase diagrams and energies for each condition
    pd_A, energy_per_atom_A, energy_above_hull_A = calculate_phase_diagram_condition_A(
        all_entries_A, entriesGases_A, locked_Chem_Potential_A, TestMat_entry_A
    )
    pd_C, energy_per_atom_C, energy_above_hull_C = calculate_phase_diagram_condition_C(
        all_entries_C, entriesGases_C, locked_Chem_Potential_C, TestMat_entry_C
    )
    pd_X, energy_per_atom_X, energy_above_hull_X = calculate_phase_diagram_condition_X(
        entriesTotal_X, entriesGases_X, locked_Chem_Potential_X, TestMat_entry_A.composition, TestMat_Ener
    )

    # Print results for condition A (Hydrogen-rich)
    print("Condition A (Hydrogen-rich):")
    print("Phase diagram:", pd_A)
    print("Energy per atom:", energy_per_atom_A)
    print("Energy above hull:", energy_above_hull_A)

    # Print results for condition C (Oxygen-rich)
    print("\nCondition C (Oxygen-rich):")
    print("Phase diagram:", pd_C)
    print("Energy per atom:", energy_per_atom_C)
    print("Energy above hull:", energy_above_hull_C)

    #Print results for condition X (CO2-rich)
    print("\nCondition X (CO2-rich):")
    print("Phase diagram:", pd_X)
    print("Energy per atom:", energy_per_atom_X)
    print("Energy above hull:", energy_above_hull_X)

if __name__ == "__main__":
    main()
