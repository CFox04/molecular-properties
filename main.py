from flask import Flask, render_template, request, abort
from rdkit.Chem import MolFromSmiles, Descriptors, Draw

app = Flask(__name__)


@app.route('/')
def home():
    return render_template('index.html')

@app.route('/smiles', methods=['POST'])
def smiles():
    smiles_str = request.form.get("smiles")
    molecule = MolFromSmiles(smiles_str)

    if molecule is None:
        return abort(400)

    svg = Draw.MolsToGridImage([molecule], useSVG=True)
    return {
        'svg': svg,
        'molProperties': {
            'Number of Radical Electrons': Descriptors.NumValenceElectrons(molecule),
            'Number of Valence Electrons': Descriptors.NumValenceElectrons(molecule),
            'Average Molecular Weight': Descriptors.MolWt(molecule),
            'Average Molecular Weight (ignoring Hydrogens)': Descriptors.HeavyAtomMolWt(molecule),
            'Exact Molecular Weight': Descriptors.ExactMolWt(molecule),
            'LogP': Descriptors.MolLogP(molecule),
            'SMILES': smiles_str
        }
    }, 200

