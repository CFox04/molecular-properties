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
    wrapper.className = 'option-item custom-control custom-checkbox mb-3'

    let checkbox = document.createElement('input');
    checkbox.type = 'checkbox';
    checkbox.checked = true;
    checkbox.name = name;
    checkbox.id = name;
    checkbox.classList.add('custom-control-input');
    
    let label = document.createElement('label');
    label.innerText = name;
    label.htmlFor = name;
    label.classList.add('custom-control-label');

    wrapper.append(checkbox);
    wrapper.append(label); 
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