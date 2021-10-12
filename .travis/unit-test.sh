#!/bin/bash

set -e
errors=0

# Run unit tests
python Magphi/Magphi_test.py || {
    echo "'python python/Magphi/Magphi_test.py' failed"
    let errors+=1
}

# Check program style
pylint -E Magphi/*.py || {
    echo 'pylint -E Magphi/*.py failed'
    let errors+=1
}

[ "$errors" -gt 0 ] && {
    echo "There were $errors errors found"
    exit 1
}

echo "Ok : Python specific tests"
