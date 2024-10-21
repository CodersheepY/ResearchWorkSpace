import os
import numpy as np
from ase.io.vasp import read_vasp_xml
from ase.dft.bandgap import bandgap

def get_bandgap(direct=False, output="bandgap_output.txt"):
    """
    Function to calculate the band gap from VASP output using ASE.

    Returns:
    gap: float
        The band gap value in eV.
    p1: tuple
        Tuple with indices of the valence band maximum.
    p2: tuple
        Tuple with indices of the conduction band minimum.
    """
    # read the vasprun.xml file
    vasp_xml_path = "vasprun.xml"
    calc = read_vasp_xml(vasp_xml_path)
    
    # calculate the band gap
    gap, p1, p2 = bandgap(calc, direct=direct)

    result = ""
    # Display and prepare the results for output
    if gap > 0:
        result = f"Band Gap: {gap:.2f} eV\nValence band max (v): {p1}\nConduction band min (c): {p2}\n"
        print(result)
    else:
        result = "The material is metallic (no band gap).\n"
        print(result)

    # Write the result to the output file
    with open(output, 'w') as f:
        f.write(result)
    
    print(f"Results written to {output}")
    
    return gap, p1, p2
if __name__ == "__main__":
    get_bandgap()
