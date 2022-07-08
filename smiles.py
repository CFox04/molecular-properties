from rdkit.Chem import MolFromSmiles, Descriptors, Draw

MOLECULE_PROPERTIES = {
    'Number of Radical Electrons': Descriptors.NumValenceElectrons,
    'Number of Valence Electrons': Descriptors.NumValenceElectrons,
    'Average Molecular Weight': Descriptors.MolWt,
    'Average Molecular Weight (ignoring hydrogens)': Descriptors.HeavyAtomMolWt,
    'Exact Molecular Weight': Descriptors.ExactMolWt,
    'LogP': Descriptors.MolLogP,
}

def get_molecule_data_from_smiles(smiles_str, options):
    molecule = MolFromSmiles(smiles_str)

    if molecule is None:
        return None

    return {
        'svg': Draw.MolsToGridImage([molecule], useSVG=True),
        'SMILES': smiles_str,
        'molProperties': {
            key:round(value(molecule), int(options['precision'])) for (key, value) in MOLECULE_PROPERTIES.items() if key in options['properties']
        }
    }
