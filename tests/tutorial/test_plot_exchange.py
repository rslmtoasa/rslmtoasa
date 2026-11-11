"""Parser/CLI checks; these do not validate the scientific exchange calculation."""
import importlib.util
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / 'tools/tutorial/plot_exchange.py'
spec = importlib.util.spec_from_file_location('plot_exchange', SCRIPT)
plot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plot)


def test_columns_and_fortran_exponents(tmp_path):
    # Catches swapped J/distance columns and unsupported Fortran exponents.
    path = tmp_path / 'jij.out'
    path.write_text('# types displacement J distance\n1 2 0.3 0.4 0 -2D-1 0.5\n')
    np.testing.assert_array_equal(plot.read_exchange(path), [[1, 2, .3, .4, 0, -.2, .5]])


@pytest.mark.parametrize('row', [
    '', '1 2 0 0 1 2', '1 2 0 0 1 2 1 3',
    '1.5 2 0 0 1 2 1', '1 0 0 0 1 2 1',
    '1 2 0 0 1 nan 1', '1 2 0 0 1 2 2',
])
def test_bad_rows(tmp_path, row):
    # Catches empty/truncated/wrong-format files, invalid types, NaNs and bad geometry.
    path = tmp_path / 'jij.out'
    path.write_text(row)
    with pytest.raises(ValueError):
        plot.read_exchange(path)


def test_filter_both_endpoints_and_preserve_directions(tmp_path):
    # Catches one-endpoint filtering and accidental removal of reverse directed rows.
    path = tmp_path / 'jij.out'
    path.write_text('1 2 .3 .4 0 -.2 .5\n2 1 -.3 -.4 0 -.2 .5\n1 3 0 0 1 9 1\n')
    figure, csv = tmp_path / 'plot.png', tmp_path / 'rows.csv'
    command = [sys.executable, str(SCRIPT), str(path), '--types', '1', '2',
               '--output', str(figure), '--csv', str(csv)]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert 'Selected 2 of 3' in result.stdout
    assert figure.read_bytes().startswith(b'\x89PNG')
    np.testing.assert_array_equal(np.loadtxt(csv, delimiter=',', skiprows=1)[:, :2], [[1, 2], [2, 1]])
    # Catches output overwrites on reruns.
    assert subprocess.run(command, capture_output=True).returncode != 0
