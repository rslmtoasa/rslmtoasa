"""
Module: test_comparison

This module contains functions to read Fortran-style namelist files and perform comparison tests
between input and reference namelists.

Usage:
    Run the tests using pytest:
        $ pytest -v test_comparison.py
"""

import pytest
from f90nml import read


def read_namelist(filename):
    """
    Read a Fortran-style namelist file and return the namelist as a dictionary.

    Parameters:
        filename (str): The path to the namelist file.

    Returns:
        dict: A dictionary containing the namelist data.
    """
    return read(filename)


def compare_namelists(input_namelist, ref_namelist, section, keys, tol=1e-6):
    """
    Compare the values in input_namelist with the reference namelist (ref_namelist).

    The function compares specified keys within a section of the input and reference namelists.
    It uses pytest.approx to allow for approximate floating-point comparisons with a specified 
    tolerance.

    Parameters:
        input_namelist (dict): The input namelist dictionary.
        ref_namelist (dict): The reference namelist dictionary.
        section (str): The section in the namelists to compare.
        keys (list): The keys to compare within the specified section.
        tol (float, optional): The tolerance value for floating-point comparisons. Default is 1e-6.

    Raises:
        AssertionError: If any of the specified keys in the input namelist do not match the
                        reference values within the given tolerance.

    Example:
        input_namelist = {'par': {'etot': 10.0, 'ws_r': 0.5, 'vmad': 2.5}}
        ref_namelist = {'par': {'etot': 10.000001, 'ws_r': 0.499999, 'vmad': 2.499999}}
        compare_namelists(input_namelist, ref_namelist, 'par', ['etot', 'ws_r', 'vmad'], tol=1e-6)
    """
    for key in keys:
        input_value = input_namelist[section][key]
        ref_value = ref_namelist[section][key]
        assert input_value == pytest.approx(ref_value, abs=tol), \
            f"Assertion failed for key: {key}. \
                Input value: {input_value}, Reference value: {ref_value}"


def test_data_comparison():
    """
    Test data comparison using pytest.

    This function reads the input and reference namelists from files 'data.nml' and 'ref.nml',
    respectively. It then performs a comparison for a specified section and keys using the
    compare_namelists function.

    To add more test cases, modify the 'section' and 'test_keys' variables as needed.

    Example:
        test_data_comparison()  # Runs the comparison test
    """
    input_namelist = read_namelist('data.nml')
    ref_namelist = read_namelist('ref.nml')
    section = 'par'
    test_keys = ['etot', 'ws_r', 'vmad']
    compare_namelists(input_namelist, ref_namelist, section, test_keys)


if __name__ == '__main__':
    """
    Run the tests using pytest.

    This code block allows you to execute the tests directly by running the script,
    invoking pytest with the '-v' flag to enable verbose output.
    """
    pytest.main(['-v', 'test_comparison.py'])
