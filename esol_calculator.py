from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, QED
from collections import namedtuple

AROMATIC_QUERY = Chem.MolFromSmarts("a")
Descriptor = namedtuple("Descriptor", "mw logp rotors ap")

def calc_ap(mol):
    """
    Calculate aromatic proportion #aromatic atoms/#atoms total
    :param mol: input molecule
    :return: aromatic proportion
    """
    matches = mol.GetSubstructMatches(AROMATIC_QUERY)
    return len(matches) / mol.GetNumAtoms()

def calc_esol_descriptors(mol):
    """
    Calcuate mw,logp,rotors and aromatic proportion (ap)
    :param mol: input molecule
    :return: named tuple with descriptor values
    """
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    rotors = Lipinski.NumRotatableBonds(mol)
    ap = calc_ap(mol)
    return Descriptor(mw=mw, logp=logp, rotors=rotors, ap=ap)

def calc_esol_orig(mol):
    """
    Original parameters from the Delaney paper, just here for comparison
    :param mol: input molecule
    :return: predicted solubility
    """
    # just here as a reference don't use this!
    intercept = 0.16
    coef = {"logp": -0.63, "mw": -0.0062, "rotors": 0.066, "ap": -0.74}
    desc = calc_esol_descriptors(mol)
    esol = intercept + coef["logp"] * desc.logp + coef["mw"] * desc.mw + coef["rotors"] * desc.rotors \
            + coef["ap"] * desc.ap
    return esol

def calc_esol(mol):
    """
    Calculate ESOL based on descriptors in the Delaney paper, coefficients refit for the RDKit using the
    routine refit_esol below
    :param mol: input molecule
    :return: predicted solubility
    """
    intercept = 0.26121066137801696
    coef = {'mw': -0.0066138847738667125, 'logp': -0.7416739523408995, 'rotors': 0.003451545565957996, 'ap': -0.42624840441316975}
    desc = calc_esol_descriptors(mol)
    esol = intercept + coef["logp"] * desc.logp + coef["mw"] * desc.mw + coef["rotors"] * desc.rotors \
            + coef["ap"] * desc.ap
    return esol