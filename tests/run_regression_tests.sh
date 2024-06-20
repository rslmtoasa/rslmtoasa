#!/bin/bash
# Script for running regression tests
# Usage: ./run_regression_tests.sh <list-of-tests>
#
# If <lists-of-tests> is not given, all defined tests are performed

# Move to the script directory
cd "$(dirname "$0")"

# Define the default list of tests
default_tests=('bccFe_lanczos' 'bccFe_block' 'bccFe_chebyshev')

# Scan for arguments
tests=()
if [ $# -eq 0 ]; then
    tests=("${default_tests[@]}")
else
    # Loop through each argument in the command-line arguments list
    for arg in "$@"; do
        if [ -d "regression/$arg" ]; then
            tests+=("$arg")
        fi
    done
fi

# Check if any valid tests were found
if [ ${#tests[@]} -eq 0 ]; then
    echo "No valid tests found."
    exit 1
fi

# Loop over tests and perform actions
echo "The following tests will be performed: ${tests[*]}"
for test in "${tests[@]}"; do
    echo "Running test: $test"
    # Uncomment the following lines to execute the test scripts
    cd "regression/$test"
    ./oneliner.sh
    cd ../..
done

