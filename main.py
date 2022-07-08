from flask import Flask, render_template, request, abort
from smiles import get_molecule_data_from_smiles, MOLECULE_PROPERTIES

app = Flask(__name__)

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/mol-properties', methods=['GET'])
def mol_properties():
    return {'properties': list(MOLECULE_PROPERTIES.keys())}, 200

@app.route('/smiles', methods=['POST'])
def smiles():
    data = get_molecule_data_from_smiles(request.json.get('smiles'), request.json.get('options'))

    if data is None:
        return abort(400)
    
    return data
