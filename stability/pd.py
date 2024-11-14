#!/usr/bin/env python

import logging
from typing import Dict, List, Tuple, Optional
import warnings

from pymatgen.core import Composition
from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram
from pymatgen.entries.computed_entries import ComputedEntry  # Correct import
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from mp_api.client import MPRester
import warnings

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def suppress_pymatgen_warnings():
    """
    抑制pymatgen关于氧化态猜测的特定警告。
    """
    warnings.filterwarnings(
        "ignore",
        message="Failed to guess oxidation states for Entry .*",
        category=UserWarning,
    )

def initialize_gas_conditions() -> Dict[str, Dict]:
    """
    初始化不同环境条件下的气体能量和化学势。
    
    Returns:
        Dict[str, Dict]: 包含每个条件下气体条目和化学势的字典。
    """
    gas_conditions = {
        "A": {  # 氢气富集环境
            "energies": {"H": -4.024, "O": -8.006},
            "composition": {"H2": "H2", "O2": "O2", "H2O": "H2O"},
            "eliminate": ["H2", "O2", "H2O"],
        },
        "C": {  # 氧气富集环境
            "energies": {"H": -4.997, "O": -6.166},
            "composition": {"H2": "H2", "O2": "O2", "H2O": "H2O"},
            "eliminate": ["H2", "O2", "H2O"],
        },
        "X": {  # 二氧化碳富集环境
            "energies": {"O": -6.166, "C": -20.232},  # 仅元素
            "composition": {"O2": "O2", "CO": "CO", "CO2": "CO2"},
            "eliminate": ["CO", "CO2", "O2"],
            "CO2_energy": -25.556,  # 单独为CO2定义能量
        },
    }
    return gas_conditions

def create_gas_entries(condition: str, gas_conditions: Dict[str, Dict]) -> Tuple[List[ComputedEntry], Dict[str, float]]:
    """
    根据条件创建气体条目和化学势。
    
    Args:
        condition (str): 环境条件 ('A', 'C', 或 'X')。
        gas_conditions (Dict[str, Dict]): 不同条件下的气体配置。
    
    Returns:
        Tuple[List[ComputedEntry], Dict[str, float]]: 气体条目和化学势。
    """
    energies = gas_conditions[condition]["energies"]
    compositions = gas_conditions[condition]["composition"]

    # 定义组合
    comp_objects = {k: Composition(v) for k, v in compositions.items()}

    gas_entries = []
    if condition in ["A", "C"]:
        # 条件A和C
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
        # 条件X
        gas_entries.append(
            ComputedEntry("O2", energies["O"] * comp_objects["O2"].num_atoms)
        )
        gas_entries.append(ComputedEntry("CO", energies["C"]))
        # 单独为CO2赋值
        co2_energy = gas_conditions[condition].get("CO2_energy", -25.556)
        gas_entries.append(ComputedEntry("CO2", co2_energy))

    # 仅元素的化学势
    chem_potentials = {}
    for element, energy in energies.items():
        chem_potentials[element] = energy

    logger.debug(f"Created gas entries for condition {condition}: {gas_entries}")
    logger.debug(f"Chemical potentials for condition {condition}: {chem_potentials}")

    return gas_entries, chem_potentials

def fetch_material_entries(api_key: str, elements: List[str], additional_elements: List[str]) -> List[ComputedEntry]:
    """
    从Materials Project获取并处理材料条目。
    
    Args:
        api_key (str): Materials Project的API密钥。
        elements (List[str]): 测试材料中的元素。
        additional_elements (List[str]): 化学系统中额外的元素。
    
    Returns:
        List[ComputedEntry]: 处理后的兼容性条目列表。
    """
    chemsys = list(set(elements + additional_elements))
    logger.info(f"Fetching entries for chemical system: {chemsys}")

    with MPRester(api_key) as mpr:
        raw_entries = mpr.get_entries_in_chemsys(chemsys)
        logger.info(f"Fetched {len(raw_entries)} raw entries from Materials Project.")

    compat = MaterialsProject2020Compatibility()
    processed_entries = compat.process_entries(raw_entries)
    logger.info(f"Processed entries: {len(processed_entries)} compatible entries.")

    return processed_entries

def calculate_phase_diagram(entries: List[ComputedEntry], chem_potentials: Dict[str, float]) -> GrandPotentialPhaseDiagram:
    """
    计算相图。
    
    Args:
        entries (List[ComputedEntry]): 条目列表。
        chem_potentials (Dict[str, float]): 化学势。
    
    Returns:
        GrandPotentialPhaseDiagram: 计算得到的相图。
    """
    try:
        pd = GrandPotentialPhaseDiagram(entries, chem_potentials)
        return pd
    except Exception as e:
        logger.error(f"Failed to create phase diagram: {e}")
        raise

def main():
    # 抑制特定的pymatgen警告
    suppress_pymatgen_warnings()

    # 配置参数
    api_key = "kzum4sPsW7GCRwtOqgDIr3zhYrfpaguK"  # 请替换为您的实际API密钥
    test_material_composition = "Ba8Fe8O24"
    test_material_energy = -268.71584216

    # 初始化气体条件
    gas_conditions = initialize_gas_conditions()

    # 定义测试材料的元素
    material_elements = [str(el) for el in Composition(test_material_composition).elements]

    # 准备每个条件下的条目和化学势
    conditions = ["A", "C", "X"]
    results = {}

    for condition in conditions:
        logger.info(f"Processing condition {condition} ({'Hydrogen-rich' if condition == 'A' else 'Oxygen-rich' if condition == 'C' else 'CO2-rich'})")

        # 创建气体条目和化学势
        gas_entries, chem_potentials = create_gas_entries(condition, gas_conditions)

        # 获取Materials Project的条目
        if condition in ["A", "C"]:
            additional_elements = ["H", "O"]
        elif condition == "X":
            additional_elements = ["C", "O"]

        material_entries = fetch_material_entries(api_key, material_elements, additional_elements)

        # 排除特定的气体物种
        eliminate = gas_conditions[condition]["eliminate"]
        filtered_entries = [e for e in material_entries if e.composition.reduced_formula not in eliminate]
        logger.debug(f"Filtered entries for condition {condition}: {len(filtered_entries)} entries.")

        # 合并气体条目
        all_entries = filtered_entries + gas_entries

        # 添加测试材料的条目
        test_entry = ComputedEntry(test_material_composition, test_material_energy)
        all_entries.append(test_entry)

        # 计算相图
        pd = calculate_phase_diagram(all_entries, chem_potentials)

        # 计算能量相关参数
        energy_per_atom = test_entry.energy / test_entry.composition.num_atoms
        gpe = next((e for e in pd.entries if e.original_entry == test_entry), None)

        if gpe:
            formation_energy = pd.get_form_energy_per_atom(gpe)
            energy_above_hull = pd.get_e_above_hull(gpe)
        else:
            logger.warning(f"Test entry not found in the phase diagram for condition {condition}.")
            formation_energy = None
            energy_above_hull = None

        # 存储结果
        results[condition] = {
            "phase_diagram": pd,
            "energy_per_atom": energy_per_atom,
            "formation_energy": formation_energy,
            "energy_above_hull": energy_above_hull,
        }

        # 输出结果
        condition_desc = "Hydrogen-rich" if condition == "A" else "Oxygen-rich" if condition == "C" else "CO2-rich"
        print(f"\nCondition {condition} ({condition_desc}):")
        print(f"Phase diagram: {pd}")
        print(f"Energy per atom: {energy_per_atom:.6f} eV")
        print(f"Formation energy: {formation_energy:.6f} eV" if formation_energy is not None else "Formation energy: N/A")
        print(f"Energy above hull: {energy_above_hull:.6f} eV" if energy_above_hull is not None else "Energy above hull: N/A")

if __name__ == "__main__":
    main()
