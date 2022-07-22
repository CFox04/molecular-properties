import { displayError } from "./error.js";

const optionsForm = document.getElementById('search-options');

export function getOptions() {
    const optionsFormData = new FormData(optionsForm);
    let options = {};

    for (const [key, value] of optionsFormData.entries()) {
        // Ticked checkboxes have a value of 'on', need to change it to true
        options[key] = value == 'on' ? true : value
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

export function displayMoleculePropertyOptions(propNames) {
    for (const propName of propNames) {
        addMoleculePropertyOption(propName);
    }
}

fetch('/mol-properties', {
    method: 'GET'
})
.then((response) => (response.ok ? response.json() : Promise.reject(response)))
.then((data) => displayMoleculePropertyOptions(data))
.catch((err) => (displayError(err)))