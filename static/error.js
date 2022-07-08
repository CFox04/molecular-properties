const errorMessage = document.getElementById('error-message');

export function displayClientError(err) {
    errorMessage.innerHTML = `<b>Unexpected Client Error</b> ${err}`;
    throw err;
}

export function displayServerError(err) {
    if (err.status == 400) {
        errorMessage.innerHTML = 'Invalid SMILES Input!'
    } else {
        errorMessage.innerHTML = `<b>Unexpected Server Error</b> (${err.status}): ${err.statusText}`;
    }
}

export function clearErrorMessage() {
    errorMessage.innerHTML = ''
}