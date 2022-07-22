import csv
from rdkit.Chem import MolFromSmiles
from io import StringIO
from smiles import MOLECULE_PROPERTIES


def get_csv_from_smiles(smiles_list, options):
    # CSV writer expects a file object, not a string. 
    # StringIO can be used to store a string as a file-like object.
    string_file = StringIO()
    writer = csv.DictWriter(string_file, fieldnames=['SMILES', *options])
    writer.writeheader()
    for smiles in smiles_list:
        molecule = MolFromSmiles(smiles)

        row = {'SMILES': smiles}

        if molecule is None:
            row['SMILES'] = f"(invalid){smiles}"

        for key, value in MOLECULE_PROPERTIES.items():
            if key in options:
                row[key] = round(value(molecule), int(options['precision']))

        writer.writerow(row)

    return string_file.getvalue()