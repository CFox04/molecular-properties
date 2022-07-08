const moleculeWrapper = document.querySelector('.molecule-wrapper');
const moleculeProperties = document.querySelector('.molecule-properties');
const moleculeSVG = document.getElementById('molecule-svg');

const FLOAT_PRECISION = 6;

function addMolecularPropertyElement(name, value) {
    let element = document.createElement('span');
    element.class = 'molecular-property';

    if (typeof value == "number" && !Number.isInteger(value)) {
        value = value.toFixed(FLOAT_PRECISION);
    }

    element.innerHTML = `<small>${name}</small> ${value}`;
    moleculeProperties.append(element);
}

export function showMoleculeWrapper() {
    if (moleculeWrapper.className.includes('hidden')) {
        moleculeWrapper.classList.remove('hidden');
    }
}
export function hideMoleculeWrapper() {
    if (!moleculeWrapper.className.includes('hidden')) {
        moleculeWrapper.classList.add('hidden');
    }
}

export function displayMoleculeCard(moleculeData) {
    showMoleculeWrapper();

    moleculeSVG.innerHTML = moleculeData.svg;

    moleculeProperties.innerHTML = '';
    for (const [propName, value] of Object.entries(moleculeData.molProperties)) {
        addMolecularPropertyElement(propName, value);
    }
}