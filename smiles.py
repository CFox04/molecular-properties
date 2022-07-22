from rdkit.Chem import MolFromSmiles, Descriptors, Draw, Crippen, Lipinski, QED
from esol_calculator import calc_ap, calc_esol, calc_esol_descriptors, calc_esol_orig

MOLECULE_PROPERTIES = {
    'Number of Radical Electrons': Descriptors.NumValenceElectrons,
    'Number of Valence Electrons': Descriptors.NumValenceElectrons,
    'Average Molecular Weight': Descriptors.MolWt,
    'Average Molecular Weight (ignoring hydrogens)': Descriptors.HeavyAtomMolWt,
    'Exact Molecular Weight': Descriptors.ExactMolWt,
    'LogP': Crippen.MolLogP,
    'Aromatic Proportion': calc_ap,
    'Number of Rotatable Bonds': Lipinski.NumRotatableBonds,
    'Water Solubility (ESOL)': calc_esol,
    'Water Solubility (ESOL Original Paper)': calc_esol_orig,
    'Polar Surface Area': lambda mol: QED.properties(mol).PSA,
}

def get_molecule_data_from_smiles(smiles_str, options):
    molecule = MolFromSmiles(smiles_str)

    if molecule is None:
        return None

    return {
        'svg': Draw.MolsToGridImage([molecule], useSVG=True, molsPerRow=1),
        'SMILES': smiles_str,
        'molProperties': {
            key:round(value(molecule), int(options['precision'])) for (key, value) in MOLECULE_PROPERTIES.items() if key in options['properties']
        }
    }
