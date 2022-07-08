import { hideMoleculeWrapper, displayMoleculeCard } from "./molecule.js";
import { clearErrorMessage, displayClientError, displayServerError } from "./error.js";

const form = document.getElementById('smiles-form')

function resetElements() {
    hideMoleculeWrapper();
    clearErrorMessage();
}

form.onsubmit = (event) => {
    const formData = new FormData(event.target);
    event.preventDefault();

    resetElements();
    
    fetch('/smiles', {
        method: 'POST',
        body: formData
    })
    .then((response) => (response.ok ? response.json() : Promise.reject(response)))
    .then((data) => displayMoleculeCard(data))
    .catch((err) => (err instanceof Response ? displayServerError(err) : displayClientError(err)))
}