import { displayError } from "./error.js";

const optionsForm = document.getElementById('search-options');

export function getOptions() {
    const optionsFormData = new FormData(optionsForm);
    let options = {'precision': optionsForm.precision.value, 'properties': []};

    for (const [key, value] of optionsFormData.entries()) {
        if (key == 'precision' || value == 'off') {
            continue;
        }
        options.properties.push(key);
    }

    return options;
}


function addMoleculePropertyOption(name) {
    let wrapper = document.createElement('div');
    wrapper.className = 'option-item'

    let checkbox = document.createElement('input');
    checkbox.type = 'checkbox';
    checkbox.checked = true;
    checkbox.name = name

    let label = document.createElement('label');
    label.innerText = name

    wrapper.append(label); 
    wrapper.append(checkbox);
    optionsForm.append(wrapper); 
}

export function displayMoleculePropertyOptions(data) {
    for (const propName of data.properties) {
        addMoleculePropertyOption(propName);
    }
}

fetch('/mol-properties', {
    method: 'GET'
})
.then((response) => (response.ok ? response.json() : Promise.reject(response)))
.then((data) => displayMoleculePropertyOptions(data))
.catch((err) => (displayError(err)))