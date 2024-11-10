#!/usr/bin/env python

import logging
from typing import Dict, List, Tuple, Optional

from pymatgen.core import Composition
from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram
from pymatgen.entries.computed_entries import ComputedEntry  # Correct import
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from mp_api.client import MPRester
import warnings

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PhaseDiagramCalculator:
    """
    A class to calculate phase diagrams under different environmental conditions
    for a specified material using pymatgen and the Materials Project API.
    """

    def __init__(
        self,
        api_key: str,
        material_composition: str,
        material_energy: float,
    ):
        """
        Initializes the PhaseDiagramCalculator with necessary parameters.

        Args:
            api_key (str): API key for the Materials Project.
            material_composition (str): Chemical formula of the test material.
            material_energy (float): Energy of the test material.
        """
        self.api_key = api_key
        self.material_composition = Composition(material_composition)
        self.material_energy = material_energy
        self.compat = MaterialsProject2020Compatibility()

        # Initialize global gas entries and chemical potentials
        self.gas_conditions = self._initialize_gas_conditions()

    @staticmethod
    def _initialize_gas_conditions() -> Dict[str, Dict]:
        """
        Initializes gas entries and chemical potentials for conditions A, C, and X.

        Returns:
            Dict[str, Dict]: A dictionary containing gas entries and chemical potentials for each condition.
        """
        gas_conditions = {
            "A": {
                "energies": {"H": -4.024, "O": -8.006},
                "composition": {"H2": "H2", "O2": "O2", "H2O": "H2O"},
                "eliminate": ["H2", "O2", "H2O"],
            },
            "C": {
                "energies": {"H": -4.997, "O": -6.166},
                "composition": {"H2": "H2", "O2": "O2", "H2O": "H2O"},
                "eliminate": ["H2", "O2", "H2O"],
            },
            "X": {
                "energies": {"O": -6.166, "C": -20.232},  # Removed 'CO2' from energies
                "composition": {"O2": "O2", "CO": "CO", "CO2": "CO2"},
                "eliminate": ["CO", "CO2", "O2"],
                "CO2_energy": -25.556,  # Added separate energy for CO2
            },
        }
        return gas_conditions

    def _create_gas_entries(
        self, condition: str
    ) -> Tuple[List[ComputedEntry], Dict[str, float]]:
        """
        Creates gas entries and locked chemical potentials for a given condition.

        Args:
            condition (str): The environmental condition ('A', 'C', or 'X').

        Returns:
            Tuple[List[ComputedEntry], Dict[str, float]]: Gas entries and chemical potentials.
        """
        energies = self.gas_conditions[condition]["energies"]
        compositions = self.gas_conditions[condition]["composition"]

        # Define compositions
        comp_objects = {k: Composition(v) for k, v in compositions.items()}

        # Create gas entries
        gas_entries = []
        if condition in ["A", "C"]:
            # For Conditions A and C
            gas_entries.append(
                ComputedEntry("H2", energies["H"] * comp_objects["H2"].num_atoms)
            )
            gas_entries.append(
                ComputedEntry("O2", energies["O"] * comp_objects["O2"].num_atoms)
            )
            gas_entries.append(
                ComputedEntry(
                    "H2O",
                    energies["H"] * comp_objects["H2"].num_atoms
                    + energies["O"] * 0.5 * comp_objects["O2"].num_atoms,
                )
            )
        elif condition == "X":
            # For Condition X
            gas_entries.append(
                ComputedEntry("O2", energies["O"] * comp_objects["O2"].num_atoms)
            )
            gas_entries.append(ComputedEntry("CO", energies["C"]))
            # Assign 'CO2' energy separately from 'energies' dict
            co2_energy = self.gas_conditions[condition].get("CO2_energy", -25.556)
            gas_entries.append(ComputedEntry("CO2", co2_energy))

        # Define chemical potentials (only elements)
        chem_potentials = {}
        for element, energy in energies.items():
            chem_potentials[element] = energy

        logger.debug(f"Created gas entries for condition {condition}: {gas_entries}")
        logger.debug(f"Chemical potentials for condition {condition}: {chem_potentials}")

        return gas_entries, chem_potentials

    def _fetch_material_entries(
        self, elements: List[str], additional_elements: List[str]
    ) -> List[ComputedEntry]:
        """
        Fetches and processes material entries from the Materials Project.

        Args:
            elements (List[str]): List of elements in the test material.
            additional_elements (List[str]): Additional elements to include in the chemical system.

        Returns:
            List[ComputedEntry]: Processed entries compatible with pymatgen's phase diagram calculations.
        """
        chemsys = list(set(elements + additional_elements))
        logger.info(f"Fetching entries for chemical system: {chemsys}")

        with MPRester(self.api_key) as mpr:
            raw_entries = mpr.get_entries_in_chemsys(chemsys)
            logger.info(f"Fetched {len(raw_entries)} raw entries from Materials Project.")

        processed_entries = self.compat.process_entries(raw_entries)
        logger.info(f"Processed entries: {len(processed_entries)} compatible entries.")

        return processed_entries

    def prepare_entries(
        self,
    ) -> Dict[str, Dict[str, List[ComputedEntry] or ComputedEntry]]:
        """
        Prepares all necessary entries for conditions A, C, and X.

        Returns:
            Dict[str, Dict[str, List[ComputedEntry] or ComputedEntry]]: Prepared entries for each condition.
        """
        entries_data = {}

        # Define additional elements based on conditions
        additional_elements = {
            "A": ["H", "O"],
            "C": ["H", "O"],
            "X": ["C", "O"],
        }

        for condition in ["A", "C", "X"]:
            gas_entries, chem_potentials = self._create_gas_entries(condition)
            eliminate = self.gas_conditions[condition]["eliminate"]

            # Fetch entries from Materials Project
            elements = [str(el) for el in self.material_composition.elements]
            entries_MP = self._fetch_material_entries(
                elements, additional_elements[condition]
            )

            # Filter out eliminated species
            filtered_entries = [
                e
                for e in entries_MP
                if e.composition.reduced_formula not in eliminate
            ]
            logger.debug(
                f"Filtered entries for condition {condition}: {len(filtered_entries)} entries."
            )

            # Add gas entries
            all_entries = filtered_entries + gas_entries

            # Create ComputedEntry for the test material
            test_material_entry = ComputedEntry(
                self.material_composition, self.material_energy
            )
            all_entries.append(test_material_entry)

            # Store in entries_data
            entries_data[condition] = {
                "all_entries": all_entries,
                "chem_potentials": chem_potentials,
                "test_material_entry": test_material_entry,
            }

        return entries_data

    def calculate_phase_diagram(
        self, condition: str, entries: List[ComputedEntry], chem_potentials: Dict[str, float], test_entry: ComputedEntry
    ) -> Tuple[GrandPotentialPhaseDiagram, float, Optional[float], Optional[float]]:
        """
        Calculates the phase diagram and related energies for a given condition.

        Args:
            condition (str): The environmental condition ('A', 'C', or 'X').
            entries (List[ComputedEntry]): List of ComputedEntry objects for the condition.
            chem_potentials (Dict[str, float]): Locked chemical potentials for the condition.
            test_entry (ComputedEntry): ComputedEntry object for the test material.

        Returns:
            Tuple[GrandPotentialPhaseDiagram, float, Optional[float], Optional[float]]: Phase diagram,
            energy per atom, formation energy, and energy above hull.
        """
        logger.info(f"Calculating phase diagram for condition {condition}.")

        # Create phase diagram
        try:
            pd = GrandPotentialPhaseDiagram(entries, chem_potentials)
        except Exception as e:
            logger.error(f"Failed to create phase diagram for condition {condition}: {e}")
            return pd, 0.0, None, None

        # Energy per atom
        energy_per_atom = test_entry.energy / test_entry.composition.num_atoms

        # Find the corresponding entry in the phase diagram
        gpe = next((e for e in pd.entries if e.original_entry == test_entry), None)

        if gpe is None:
            logger.warning(f"Test entry not found in the phase diagram for condition {condition}.")
            formation_energy = None
            energy_above_hull = None
        else:
            # Formation energy
            formation_energy = pd.get_form_energy_per_atom(gpe)

            # Energy above hull
            energy_above_hull = pd.get_e_above_hull(gpe)

        logger.debug(
            f"Condition {condition}: Energy per atom = {energy_per_atom}, "
            f"Formation energy = {formation_energy}, Energy above hull = {energy_above_hull}"
        )

        return pd, energy_per_atom, formation_energy, energy_above_hull

    def run_calculations(self):
        """
        Runs the phase diagram calculations for all conditions and prints the results.
        """
        # Prepare all entries
        entries_data = self.prepare_entries()

        results = {}

        # Calculate phase diagrams for each condition
        for condition in ["A", "C", "X"]:
            pd, energy_per_atom, formation_energy, energy_above_hull = self.calculate_phase_diagram(
                condition,
                entries_data[condition]["all_entries"],
                entries_data[condition]["chem_potentials"],
                entries_data[condition]["test_material_entry"],
            )

            results[condition] = {
                "phase_diagram": pd,
                "energy_per_atom": energy_per_atom,
                "formation_energy": formation_energy,
                "energy_above_hull": energy_above_hull,
            }

            # Print results
            print(f"\nCondition {condition} ({self._condition_description(condition)}):")
            print(f"Phase diagram: {pd}")
            print(f"Energy per atom: {energy_per_atom:.6f} eV")
            print(f"Formation energy: {formation_energy:.6f} eV" if formation_energy is not None else "Formation energy: N/A")
            print(f"Energy above hull: {energy_above_hull:.6f} eV" if energy_above_hull is not None else "Energy above hull: N/A")

    @staticmethod
    def _condition_description(condition: str) -> str:
        """
        Provides a textual description for each condition.

        Args:
            condition (str): The environmental condition ('A', 'C', or 'X').

        Returns:
            str: Description of the condition.
        """
        descriptions = {
            "A": "Hydrogen-rich",
            "C": "Oxygen-rich",
            "X": "CO2-rich",
        }
        return descriptions.get(condition, "Unknown")


def suppress_pymatgen_warnings():
    """
    Suppress specific pymatgen warnings related to oxidation state guessing.
    """
    warnings.filterwarnings(
        "ignore",
        message="Failed to guess oxidation states for Entry .*",
        category=UserWarning,
    )


def main():
    # Suppress specific pymatgen warnings
    suppress_pymatgen_warnings()

    # API key and test material parameters
    api_key = "kzum4sPsW7GCRwtOqgDIr3zhYrfpaguK"  # Replace with your actual API key
    test_material_composition = "Ba8Zr8O24"
    test_material_energy = -338.71584216

    # Initialize the calculator
    calculator = PhaseDiagramCalculator(
        api_key=api_key,
        material_composition=test_material_composition,
        material_energy=test_material_energy,
    )

    # Run the phase diagram calculations
    calculator.run_calculations()


if __name__ == "__main__":
    main()
