from rdkit import Chem
from rdkit.Chem import AllChem
import json
import sys

def mol2_to_ccd_cif(mol2_file: str, component_id: str, template_smiles: str = None):
    # Load the mol2 file as an RDKit Mol object
    mol = Chem.MolFromMol2File(mol2_file, removeHs=False)
    if mol is None:
        raise ValueError("Failed to load mol2 file or invalid mol2 format.")

    # Generate 3D coordinates if missing
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)

    # Attempt to sanitize the molecule for bond orders and aromaticity
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        print(f"Sanitization failed: {e}")

    # Optionally use a template to assign bond orders more accurately
    if template_smiles:
        template_mol = Chem.MolFromSmiles(template_smiles)
        if template_mol:
            AllChem.Compute2DCoords(template_mol)
            try:
                mol = Chem.AssignBondOrdersFromTemplate(template_mol, mol)
            except Exception as e:
                print(f"Template matching failed: {e}")
        else:
            raise ValueError("Invalid SMILES for the template molecule.")

    # Extract atom names from MOL2
    with open(mol2_file, 'r') as f:
        lines = f.readlines()

    atom_name_map = {}
    in_atom_section = False
    for line in lines:
        if line.strip().startswith("@<TRIPOS>ATOM"):
            in_atom_section = True
            continue
        if line.strip().startswith("@<TRIPOS>BOND"):
            break
        if in_atom_section:
            parts = line.split()
            if len(parts) >= 2:
                idx = int(parts[0]) - 1  # MOL2 is 1-indexed
                name = parts[1]
                atom_name_map[idx] = name

    # Convert to CCD-like mmCIF with hydrogens excluded
    cif_content = []
    cif_content.append("data_" + component_id + "\n#")
    cif_content.append("_chem_comp.id " + component_id)
    cif_content.append("_chem_comp.name ?")
    cif_content.append("_chem_comp.type lipid")
    cif_content.append("_chem_comp.formula ?")
    cif_content.append("_chem_comp.mon_nstd_parent_comp_id ?")
    cif_content.append("_chem_comp.pdbx_synonyms ?")
    cif_content.append("_chem_comp.formula_weight 747.07")
    cif_content.append("#")
    cif_content.append("loop_")
    cif_content.append("_chem_comp_atom.comp_id")
    cif_content.append("_chem_comp_atom.atom_id")
    cif_content.append("_chem_comp_atom.type_symbol")
    cif_content.append("_chem_comp_atom.charge")
    cif_content.append("_chem_comp_atom.pdbx_leaving_atom_flag")
    cif_content.append("_chem_comp_atom.pdbx_model_Cartn_x_ideal")
    cif_content.append("_chem_comp_atom.pdbx_model_Cartn_y_ideal")
    cif_content.append("_chem_comp_atom.pdbx_model_Cartn_z_ideal")

    atom_dict = {}

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            continue
        idx = atom.GetIdx()
        element = atom.GetSymbol()
        charge = 0
        coords = mol.GetConformer().GetAtomPosition(idx)
        atom_name = atom_name_map.get(idx, f"{element}{idx+1}")
        atom_dict[idx] = atom_name
        cif_content.append(f"{component_id} {atom_name} {element} {charge} N {coords.x:.3f} {coords.y:.3f} {coords.z:.3f}")

    cif_content.append("#")
    cif_content.append("loop_")
    cif_content.append("_chem_comp_bond.atom_id_1")
    cif_content.append("_chem_comp_bond.atom_id_2")
    cif_content.append("_chem_comp_bond.value_order")
    cif_content.append("_chem_comp_bond.pdbx_aromatic_flag")

    bonds_to_add = []
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtomIdx()
        atom2 = bond.GetEndAtomIdx()
        if atom1 not in atom_dict or atom2 not in atom_dict:
            continue
        bond_order = 'SING' if bond.GetBondType() == Chem.BondType.SINGLE else 'DOUB'
        atom1_name = atom_dict[atom1]
        atom2_name = atom_dict[atom2]
        bonds_to_add.append((atom1_name, atom2_name, bond_order))

    sorted_bonds = sorted(bonds_to_add, key=lambda bond: bond[0])

    for bond in sorted_bonds:
        cif_content.append(f"{bond[0]} {bond[1]} {bond[2]} N")

    output_cif = f"{component_id}.cif"
    with open(output_cif, 'w') as cif_file:
        cif_file.write("\n".join(cif_content))

    json_content = {
        "userCCD": "\n".join(cif_content),
    }
    output_json = f"{component_id}.json"
    with open(output_json, 'w') as json_file:
        json.dump(json_content, json_file, indent=4)

    print(f"Generated mmCIF: {output_cif}")
    print(f"Generated JSON: {output_json}")

if __name__ == "__main__":
    mol2_file = sys.argv[1]
    component_id = sys.argv[2]
    template_smiles = None
    mol2_to_ccd_cif(mol2_file, component_id, template_smiles)

