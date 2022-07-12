const errorMessage = document.getElementById('error-message');

function displayClientError(err) {
    errorMessage.innerHTML = `<b>Unexpected Client Error</b> ${err}`;
    throw err;
}

function displayServerError(err) {
    if (err.status == 400) {
        errorMessage.innerHTML = 'Invalid SMILES Input!';
    } else {
        errorMessage.innerHTML = `<b>Unexpected Server Error</b> (${err.status}): ${err.statusText}`;
    }
}

export function clearErrorMessage() {
    errorMessage.innerHTML = '';
}

export function displayError(err) {
    if (err instanceof Response) { 
        displayServerError(err) 
    } else {
        displayClientError(err)
    } 
}