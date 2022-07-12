import { getOptions } from '/static/options.js';
import { displayError, clearErrorMessage } from "./error.js";

const form = document.getElementById('multi-smiles-form')
const textArea = document.getElementById('multi-smiles')

function downloadCSV(blob) {
    let url = window.URL.createObjectURL(blob);
    window.location.assign(url);
    URL.revokeObjectURL(url);
}

form.onsubmit = (event) => {
    event.preventDefault();

    clearErrorMessage();
    
    fetch('/smiles-csv', {
        method: 'POST',
        body: JSON.stringify({'smiles': textArea.value.replace(/\s/g, '').split(','), 'options': getOptions()}),
        headers: {
            'Content-type': 'application/json; charset=UTF-8',
        }
    })
    .then(async (response) => (response.ok ? response.blob() : Promise.reject(response)))
    .then((blob) => downloadCSV(blob))
    .catch((err) => (displayError(err)))
}