Flask app that gets various properties of a molecule provided in SMILES format.

## Running the program

1. Install [anaconda](https://www.anaconda.com/products/distribution)
2. Create conda environment from the provided .yaml ``conda env create --file environment.yaml``
3. Activate the environment ``conda activate cenv``
4. Set the FLASK_APP [environment variable](https://flask.palletsprojects.com/en/2.1.x/config/#environment-and-debug-features) (Powershell): ``$env:FLASK_ENV = "development"``
5. Run the flask app ``flask run``