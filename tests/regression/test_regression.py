"""
Module: my_data_processing

This module contains functions to read Fortran-style namelist files, process the data, and export it to YAML format.
It also includes a test function to check data using num_regression.

Usage:
    To read and parse data from a namelist file:
        data = parse_inputdata('Fe_out.nml')

    To export data to a YAML file:
        data = {'par': {'key1': value1, 'key2': value2, ...}}
        export_data(data, write_to_file=True)

    To perform data comparison tests:
        pytest -v test_data_processing.py
"""
from collections import OrderedDict
import yaml
import numpy as np
from f90nml import read


def read_namelist_to_dict(filename):
    """
    Read a Fortran-style namelist file and return the namelist as a dictionary.

    Parameters:
        filename (str): The path to the namelist file.

    Returns:
        dict: A dictionary containing the namelist data.
    """
    nml = read(filename)
    odict = nml.todict()

    try:
        del odict['par']['_start_index']
    except KeyError:
        pass

    yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping(
        'tag:yaml.org,2002:map', data.items()))

    outdict = yaml.load(yaml.dump(odict), Loader=yaml.Loader)

    return outdict


def export_data(indict, write_to_file=False):
    """
    Export data in the form of a dictionary to a YAML file.

    Parameters:
        indict (dict): The input dictionary to be exported.
        write_to_file (bool, optional): Whether to write the YAML content to a file.
                                       If True, the content is written to 'out.yaml'. Default is False.

    Example:
        data = {'par': {'key1': value1, 'key2': value2, ...}}
        export_data(data, write_to_file=True)  # Exports data to 'out.yaml'
    """
    with open('out.yaml', 'w') as yfile:
        yaml.dump(indict, yfile)

    yaml_string = yaml.dump(indict)

    print(yaml_string)


def parse_inputdata(filename):
    """
    Parse input data from a namelist file and return a processed dictionary.

    Parameters:
        filename (str): The path to the namelist file.

    Returns:
        dict: A dictionary containing the processed namelist data.

    Example:
        data = parse_inputdata('Fe_out.nml')
    """
    i_dict = read_namelist_to_dict(filename)['par']

    for k, obj in i_dict.items():
        if not isinstance(obj, np.ndarray):
            arr = np.atleast_1d(np.asarray(obj))
            if np.issubdtype(arr.dtype, np.number):
                i_dict[k] = arr

    for key in i_dict.keys():
        tarray = np.array(i_dict[key]).ravel()
        tarray = tarray[tarray != None]
        i_dict[key] = float(tarray[-1])

    return i_dict


def test_data(num_regression):
    """
    Test the data using num_regression.

    Parameters:
        num_regression: The num_regression fixture provided by pytest-regressions.

    Example:
        test_data(num_regression)
    """
    filename = 'data.nml'
    data = parse_inputdata(filename)
    num_regression.check(data, default_tolerance=dict(atol=1e-6, rtol=1e-6))


if __name__ == '__main__':
    FILENAME = 'data.nml'
    data = parse_inputdata(FILENAME)
    export_data(data)
