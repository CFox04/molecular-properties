from typing import OrderedDict
from rdkit.Chem import MolFromSmiles, Descriptors, Draw, Crippen, Lipinski, QED
from esol_calculator import calc_ap, calc_esol, calc_esol_descriptors, calc_esol_orig

MOLECULE_PROPERTIES = OrderedDict({
    'Number of Radical Electrons': Descriptors.NumValenceElectrons,
    'Number of Valence Electrons': Descriptors.NumValenceElectrons,
    'Average Molecular Weight': Descriptors.MolWt,
    'Average Molecular Weight (ignoring hydrogens)': Descriptors.HeavyAtomMolWt,
    'Exact Molecular Weight': Descriptors.ExactMolWt,
    'LogP': Crippen.MolLogP,
    'Aromatic Proportion': calc_ap,
    'Number of Rotatable Bonds': Lipinski.NumRotatableBonds,
    'Water Solubility (ESOL Original Paper)': calc_esol_orig,
    'Water Solubility (ESOL)': calc_esol,
    'Polar Surface Area': lambda mol: QED.properties(mol).PSA,
})

def get_molecule_data_from_smiles(smiles_str, options):
    molecule = MolFromSmiles(smiles_str)

    if molecule is None:
        return None
    
    def get_value(value_func):
        return round(value_func(molecule), int(options['precision']))

    return {
        'svg': Draw.MolsToGridImage([molecule], useSVG=True, molsPerRow=1),
        'SMILES': smiles_str,
        # Use list of tuples instead of dict to maintain key order
        'molProperties': [(prop_name, get_value(value_func)) for (prop_name, value_func) in MOLECULE_PROPERTIES.items() if prop_name in options]
    }
